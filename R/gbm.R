.onAttach <- function(lib, pkg)
{
   vers <- library(help=gbm)$info[[1]]
   vers <- vers[grep("Version:",vers)]
   vers <- rev(strsplit(vers," ")[[1]])[1]
   packageStartupMessage(paste("Loaded gbm",vers))
}

gbm <- function(formula = formula(data),
                distribution = "bernoulli",
                data = list(),
                weights,
                var.monotone = NULL,
                n.trees = 100,
                interaction.depth = 1,
                n.minobsinnode = 10,
                shrinkage = 0.001,
                bag.fraction = 0.5,
                train.fraction = 1.0,
                cv.folds=0,
                keep.data = TRUE,
                verbose = TRUE,
                class.stratify.cv)
{
   theCall <- match.call()

   mf <- match.call(expand.dots = FALSE)
   m <- match(c("formula", "data", "weights", "offset"), names(mf), 0)
   mf <- mf[c(1, m)]
   mf$drop.unused.levels <- TRUE
   mf$na.action <- na.pass
   mf[[1]] <- as.name("model.frame")
   m <- mf
   mf <- eval(mf, parent.frame())
   Terms <- attr(mf, "terms")

   y <- model.response( mf )
   # If distribution is not given, try to guess it
   if ( missing( distribution ) ){
      if ( length( unique( y ) ) == 2 ){ distribution <- "bernoulli" }
      else if ( class( y ) == "Surv" ){ distribution <- "coxph" }
      else if ( is.factor( y ) ){ distribution <- "multinomial" }
      else{
         distribution <- "gaussian"
      }
      cat( paste( "Distribution not specified, assuming", distribution, "...\n" ) )
   }

   #   if ( length( distribution ) == 1 && distribution != "multinomial" ){
   #      y <- model.response(mf, "numeric")
   #   }
   #   else { y <- model.response( mf ) }

   w <- model.weights(mf)
   offset <- model.offset(mf)

   var.names <- attributes(Terms)$term.labels
   x <- model.frame(terms(reformulate(var.names)),
                    data,
                    na.action=na.pass)

   # get the character name of the response variable
   response.name <- as.character(formula[[2]])
   if(is.character(distribution)) distribution <- list(name=distribution)

   if ( missing( class.stratify.cv ) ){
      if ( distribution$name == "multinomial" ){ class.stratify.cv <- TRUE }
      else class.stratify.cv <- FALSE
   }
   else {
      if ( !is.element( distribution$name, c( "bernoulli", "multinomial" ) ) ){
         warning("You can only use class.stratify.cv when distribution is bernoulli or multinomial. Ignored.")
         class.stratify.cv <- FALSE
      }
   }

   # groups (for pairwise distribution only)
   group      <- NULL
   num.groups <- 0

   # determine number of training instances
   if(distribution$name=="coxph")
   {
      nTrain <- floor(train.fraction*nrow(y))
   }
   else if (distribution$name != "pairwise")
   {
      nTrain <- floor(train.fraction*length(y))
   }
   else
   {
      # distribution$name == "pairwise":
      # Sampling is by group, so we need to calculate them here
      distribution.group <- distribution[["group"]]
      if (is.null(distribution.group))
      {
         stop("For pairwise regression, the distribution parameter must be a list with a parameter 'group' for the a list of the column names indicating groups, for example list(name=\"pairwise\",group=c(\"date\",\"session\",\"category\",\"keywords\")).")
      }

      # Check if group names are valid
      i <- match(distribution.group, colnames(data))
      if (any(is.na(i)))
      {
         stop("Group column does not occur in data: ", distribution.group[is.na(i)])
      }

      # Construct group index
      group <- factor(do.call(paste, c(data[,distribution.group, drop=FALSE], sep=":")))

      # Check that weights are constant across groups
      if ((!missing(weights)) && (!is.null(weights)))
      {
         w.min <- tapply(w, INDEX=group, FUN=min)
         w.max <- tapply(w, INDEX=group, FUN=max)

         if (any(w.min != w.max))
         {
            stop("For distribution 'pairwise', all instances for the same group must have the same weight")
         }

         # Normalize across groups
         w <- w * length(w.min) / sum(w.min)
      }

      # Shuffle groups, to remove bias when splitting into train/test set and/or CV folds
      perm.levels  <- levels(group)[sample(1:nlevels(group))]
      group        <- factor(group, levels=perm.levels)

      # The C function expects instances to be sorted by group and descending by target
      ord.group    <- order(group, -y)
      group        <- group[ord.group]
      y            <- y[ord.group]
      x            <- x[ord.group,,drop=FALSE]
      w            <- w[ord.group]

      # Split into train and validation set, at group boundary
      num.groups.train <- max(1, round(train.fraction * nlevels(group)))

      # include all groups up to the num.groups.train
      nTrain           <- max(which(group==levels(group)[num.groups.train]))
      Misc             <- group
   } # close if(distribution$name=="coxph") ...

   cv.error <- NULL
   if(cv.folds>1)
   {
      i.train <- 1:nTrain

      if ( distribution$name %in% c( "bernoulli", "multinomial" ) & class.stratify.cv )
      {
         nc <- table(y[i.train]) # Number in each class
         uc <- names(nc)
         if ( min( nc ) < cv.folds ){
            stop( paste("The smallest class has only", min(nc), "objects in the training set. Can't do", cv.folds, "fold cross-validation."))
         }
         cv.group <- vector( length = length( i.train ) )
         for ( i in 1:length( uc ) ){
            cv.group[ y[i.train] == uc[i] ] <- sample( rep( 1:cv.folds , length = nc[i] ) )
         }
      }
      else if (distribution$name == "pairwise")
      {
         # Split into CV folds at group boundaries
         s <- sample(rep(1:cv.folds, length=nlevels(group)))
         cv.group <- s[as.integer(group[i.train])]
      }
      else
      {
         cv.group <- sample(rep(1:cv.folds, length=length(i.train)))
      }

      cv.error <- rep(0, n.trees)
      for(i.cv in 1:cv.folds)
      {
         if(verbose) cat("CV:",i.cv,"\n")

         i <- order(cv.group==i.cv)

         gbm.obj <- gbm.fit(x[i.train,,drop=FALSE][i,,drop=FALSE],
                            y[i.train][i],
                            offset = offset[i.train][i],
                            distribution = distribution,
                            w = if(is.null(w)) logical(0) else w[i.train][i],
                            var.monotone = var.monotone,
                            n.trees = n.trees,
                            interaction.depth = interaction.depth,
                            n.minobsinnode = n.minobsinnode,
                            shrinkage = shrinkage,
                            bag.fraction = bag.fraction,
                            nTrain = length(which(cv.group!=i.cv)),
                            keep.data = FALSE,
                            verbose = verbose,
                            var.names = var.names,
                            response.name = response.name,
                            group = group[i.train][i])
         cv.error <- cv.error + gbm.obj$valid.error*sum(cv.group==i.cv)
      }
      cv.error <- cv.error/length(i.train)
   }

   gbm.obj <- gbm.fit(x,y,
                      offset = offset,
                      distribution = distribution,
                      w = w,
                      var.monotone = var.monotone,
                      n.trees = n.trees,
                      interaction.depth = interaction.depth,
                      n.minobsinnode = n.minobsinnode,
                      shrinkage = shrinkage,
                      bag.fraction = bag.fraction,
                      nTrain = nTrain,
                      keep.data = keep.data,
                      verbose = verbose,
                      var.names = var.names,
                      response.name = response.name,
                      group = group)

   gbm.obj$train.fraction <- train.fraction
   gbm.obj$Terms <- Terms
   gbm.obj$cv.error <- cv.error
   gbm.obj$cv.folds <- cv.folds
   gbm.obj$call <- theCall
   gbm.obj$m <- m

   if (distribution$name == "pairwise")
   {
      # Data has been reordered according to queries.
      # We need to permute the fitted values to correspond
      # to the original order.
      gbm.obj$ord.group <- ord.group
      gbm.obj$fit <- gbm.obj$fit[order(ord.group)]
   }

   return(gbm.obj)
}
