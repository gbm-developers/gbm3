gbm.fit <- function(x,y,
                    offset = NULL,
                    misc = NULL,
                    distribution = "bernoulli",
                    w = NULL,
                    var.monotone = NULL,
                    n.trees = 100,
                    interaction.depth = 1,
                    n.minobsinnode = 10,
                    shrinkage = 0.001,
                    bag.fraction = 0.5,
                    nTrain = NULL,
                    train.fraction = NULL,
                    keep.data = TRUE,
                    verbose = TRUE,
                    var.names = NULL,
                    response.name = "y",
                    group = NULL)
{

   if(is.character(distribution)) { distribution <- list(name=distribution) }

   cRows <- nrow(x)
   cCols <- ncol(x)

   if(nrow(x) != ifelse(class(y)=="Surv", nrow(y), length(y))) {
      stop("The number of rows in x does not equal the length of y.")
   }

   # the preferred way to specify the number of training instances is via parameter 'nTrain'.
   # parameter 'train.fraction' is only maintained for backward compatibility.

   if(!is.null(nTrain) && !is.null(train.fraction)) {
      stop("Parameters 'nTrain' and 'train.fraction' cannot both be specified")
   }
   else if(!is.null(train.fraction)) {
      warning("Parameter 'train.fraction' of gbm.fit is deprecated, please specify 'nTrain' instead")
      nTrain <- floor(train.fraction*cRows)
   }
   else if(is.null(nTrain)) {
     # both undefined, use all training data
     nTrain <- cRows
   }

   if (is.null(train.fraction)){
      train.fraction <- nTrain / cRows
   }

   if(is.null(var.names)) {
       var.names <- getVarNames(x)
   }

#   if(is.null(response.name)) { response.name <- "y" }

   # check dataset size
   if(nTrain * bag.fraction <= 2*n.minobsinnode+1) {
      stop("The dataset size is too small or subsampling rate is too large: nTrain*bag.fraction <= n.minobsinnode")
   }

   if (distribution$name != "pairwise") {
      w <- w*length(w)/sum(w) # normalize to N
   }

   # Do sanity checks
   ch <- checkMissing(x, y)
   interaction.depth <- checkID(interaction.depth)
   w <- checkWeights(w, length(y))
   offset <- checkOffset(offset, y)

   Misc <- NA

   # setup variable types
   var.type <- rep(0,cCols)
   var.levels <- vector("list",cCols)
   for(i in 1:length(var.type))
   {
      if(all(is.na(x[,i])))
      {
         stop("variable ",i,": ",var.names[i]," has only missing values.")
      }
      if(is.ordered(x[,i]))
      {
         var.levels[[i]] <- levels(x[,i])
         x[,i] <- as.numeric(x[,i])-1
         var.type[i] <- 0
      }
      else if(is.factor(x[,i]))
      {
         if(length(levels(x[,i]))>1024)
            stop("gbm does not currently handle categorical variables with more than 1024 levels. Variable ",i,": ",var.names[i]," has ",length(levels(x[,i]))," levels.")
         var.levels[[i]] <- levels(x[,i])
         x[,i] <- as.numeric(x[,i])-1
         var.type[i] <- max(x[,i],na.rm=TRUE)+1
      }
      else if(is.numeric(x[,i]))
      {
         var.levels[[i]] <- quantile(x[,i],prob=(0:10)/10,na.rm=TRUE)
      }
      else
      {
         stop("variable ",i,": ",var.names[i]," is not of type numeric, ordered, or factor.")
      }

      # check for some variation in each variable
      if(length(unique(var.levels[[i]])) == 1)
      {
         warning("variable ",i,": ",var.names[i]," has no variation.")
      }
   }

   nClass <- 1

   if(!("name" %in% names(distribution))) {
      stop("The distribution is missing a 'name' component, for example list(name=\"gaussian\")")
   }
   supported.distributions <-
   c("bernoulli","gaussian","poisson","adaboost","laplace","coxph","quantile",
     "tdist", "multinomial", "huberized", "pairwise")

   distribution.call.name <- distribution$name

   # check potential problems with the distributions
   if(!is.element(distribution$name,supported.distributions))
   {
      stop("Distribution ",distribution$name," is not supported")
   }
   if((distribution$name == "bernoulli") && !all(is.element(y,0:1)))
   {
      stop("Bernoulli requires the response to be in {0,1}")
   }
   if((distribution$name == "huberized") && !all(is.element(y,0:1)))
   {
      stop("Huberized square hinged loss requires the response to be in {0,1}")
   }
   if((distribution$name == "poisson") && any(y<0))
   {
      stop("Poisson requires the response to be positive")
   }
   if((distribution$name == "poisson") && any(y != trunc(y)))
   {
      stop("Poisson requires the response to be a positive integer")
   }
   if((distribution$name == "adaboost") && !all(is.element(y,0:1)))
   {
      stop("This version of AdaBoost requires the response to be in {0,1}")
   }
   if(distribution$name == "quantile")
   {
      if(length(unique(w)) > 1)
      {
         stop("This version of gbm for the quantile regression lacks a weighted quantile. For now the weights must be constant.")
      }
      if(is.null(distribution$alpha))
      {
         stop("For quantile regression, the distribution parameter must be a list with a parameter 'alpha' indicating the quantile, for example list(name=\"quantile\",alpha=0.95).")
      } else
      if((distribution$alpha<0) || (distribution$alpha>1))
      {
         stop("alpha must be between 0 and 1.")
      }
      Misc <- c(alpha=distribution$alpha)
   }
   if(distribution$name == "coxph")
   {
      if(class(y)!="Surv")
      {
         stop("Outcome must be a survival object Surv(time,failure)")
      }
      if(attr(y,"type")!="right")
      {
         stop("gbm() currently only handles right censored observations")
      }
      Misc <- y[,2]
      y <- y[,1]

      # reverse sort the failure times to compute risk sets on the fly
      i.train <- order(-y[1:nTrain])
      n.test <- cRows - nTrain
      if(n.test > 0)
      {
         i.test <- order(-y[(nTrain+1):cRows]) + nTrain
      }
      else
      {
         i.test <- NULL
      }
      i.timeorder <- c(i.train,i.test)

      y <- y[i.timeorder]
      Misc <- Misc[i.timeorder]
      x <- x[i.timeorder,,drop=FALSE]
      w <- w[i.timeorder]
      if(!is.na(offset)) offset <- offset[i.timeorder]
   }
   if(distribution$name == "tdist")
   {
      if (is.null(distribution$df) || !is.numeric(distribution$df)){
         Misc <- 4
      }
      else {
         Misc <- distribution$df[1]
      }
   }
   if (distribution$name == "multinomial")
   {
      ## Ensure that the training set contains all classes
      classes <- attr(factor(y), "levels")
      nClass <- length(classes)

      if (nClass > nTrain){
         stop(paste("Number of classes (", nClass,
                    ") must be less than the size of the training set (", nTrain, ")",
                    sep = ""))
      }

      #    f <- function(a,x){
      #       min((1:length(x))[x==a])
      #    }

      new.idx <- as.vector(sapply(classes, function(a,x){ min((1:length(x))[x==a]) }, y))

      all.idx <- 1:length(y)
      new.idx <- c(new.idx, all.idx[!(all.idx %in% new.idx)])

      y <- y[new.idx]
      x <- x[new.idx, ]
      w <- w[new.idx]
      if (!is.null(offset)){
         offset <- offset[new.idx]
      }

      ## Get the factors
      y <- as.numeric(as.vector(outer(y, classes, "==")))

      ## Fill out the weight and offset
      w <- rep(w, nClass)
      if (!is.null(offset)){
         offset <- rep(offset, nClass)
      }
   } # close if (dist... == "multinomial"

   if(distribution$name == "pairwise")
   {
      distribution.metric <- distribution[["metric"]]
      if (!is.null(distribution.metric))
      {
         distribution.metric <- tolower(distribution.metric)
         supported.metrics <- c("conc", "ndcg", "map", "mrr")
         if (!is.element(distribution.metric, supported.metrics))
         {
            stop("Metric '", distribution.metric, "' is not supported, use either 'conc', 'ndcg', 'map', or 'mrr'")
         }
         metric <- distribution.metric
      }
      else
      {
         warning("No metric specified, using 'ndcg'")
         metric <- "ndcg" # default
         distribution[["metric"]] <- metric
      }

      if (any(y<0))
      {
         stop("targets for 'pairwise' should be non-negative")
      }

      if (is.element(metric, c("mrr", "map")) && (!all(is.element(y, 0:1))))
      {
         stop("Metrics 'map' and 'mrr' require the response to be in {0,1}")
      }

      # Cut-off rank for metrics
      # Default of 0 means no cutoff

      max.rank <- 0
      if (!is.null(distribution[["max.rank"]]) && distribution[["max.rank"]] > 0)
      {
         if (is.element(metric, c("ndcg", "mrr")))
         {
            max.rank <- distribution[["max.rank"]]
         }
         else
         {
            stop("Parameter 'max.rank' cannot be specified for metric '", distribution.metric, "', only supported for 'ndcg' and 'mrr'")
         }
      }

      # We pass the cut-off rank to the C function as the last element in the Misc vector
      Misc <- c(group, max.rank)

      distribution.call.name <- sprintf("pairwise_%s", metric)
   } # close if (dist... == "pairwise"

   # create index upfront... subtract one for 0 based order
   x.order <- apply(x[1:nTrain,,drop=FALSE],2,order,na.last=FALSE)-1

   x <- as.vector(data.matrix(x))
   predF <- rep(0,length(y))
   train.error <- rep(0,n.trees)
   valid.error <- rep(0,n.trees)
   oobag.improve <- rep(0,n.trees)

   if(is.null(var.monotone)) var.monotone <- rep(0,cCols)
   else if(length(var.monotone)!=cCols)
   {
      stop("Length of var.monotone != number of predictors")
   }
   else if(!all(is.element(var.monotone,-1:1)))
   {
      stop("var.monotone must be -1, 0, or 1")
   }
   fError <- FALSE

   gbm.obj <- .Call("gbm",
                    Y=as.double(y),
                    Offset=as.double(offset),
                    X=as.double(x),
                    X.order=as.integer(x.order),
                    weights=as.double(w),
                    Misc=as.double(Misc),
                    cRows=as.integer(cRows),
                    cCols=as.integer(cCols),
                    var.type=as.integer(var.type),
                    var.monotone=as.integer(var.monotone),
                    distribution=as.character(distribution.call.name),
                    n.trees=as.integer(n.trees),
                    interaction.depth=as.integer(interaction.depth),
                    n.minobsinnode=as.integer(n.minobsinnode),
                    n.classes = as.integer(nClass),
                    shrinkage=as.double(shrinkage),
                    bag.fraction=as.double(bag.fraction),
                    nTrain=as.integer(nTrain),
                    fit.old=as.double(NA),
                    n.cat.splits.old=as.integer(0),
                    n.trees.old=as.integer(0),
                    verbose=as.integer(verbose),
                    PACKAGE = "gbm")

   names(gbm.obj) <- c("initF","fit","train.error","valid.error",
                       "oobag.improve","trees","c.splits")

   gbm.obj$bag.fraction <- bag.fraction
   gbm.obj$distribution <- distribution
   gbm.obj$interaction.depth <- interaction.depth
   gbm.obj$n.minobsinnode <- n.minobsinnode
   gbm.obj$num.classes <- nClass
   gbm.obj$n.trees <- length(gbm.obj$trees) / nClass
   gbm.obj$nTrain <- nTrain
   gbm.obj$train.fraction <- train.fraction
   gbm.obj$response.name <- response.name
   gbm.obj$shrinkage <- shrinkage
   gbm.obj$var.levels <- var.levels
   gbm.obj$var.monotone <- var.monotone
   gbm.obj$var.names <- var.names
   gbm.obj$var.type <- var.type
   gbm.obj$verbose <- verbose
   gbm.obj$Terms <- NULL

   if(distribution$name == "coxph")
   {
      gbm.obj$fit[i.timeorder] <- gbm.obj$fit
   }
   ## If K-Classification is used then split the fit and tree components
   if (distribution$name == "multinomial"){
      gbm.obj$fit <- matrix(gbm.obj$fit, ncol = nClass)
      dimnames(gbm.obj$fit)[[2]] <- classes
      gbm.obj$classes <- classes

      ## Also get the class estimators
      exp.f <- exp(gbm.obj$fit)
      denom <- matrix(rep(rowSums(exp.f), nClass), ncol = nClass)
      gbm.obj$estimator <- exp.f/denom
   }

   if(keep.data)
   {
      if(distribution$name == "coxph")
      {
         # put the observations back in order
         gbm.obj$data <- list(y=y,x=x,x.order=x.order,offset=offset,Misc=Misc,w=w,
                              i.timeorder=i.timeorder)
      }
      else if ( distribution$name == "multinomial" ){
         # Restore original order of the data
         new.idx <- order( new.idx )
         gbm.obj$data <- list( y=as.vector(matrix(y, ncol=length(classes),byrow=FALSE)[new.idx,]),
                              x=as.vector(matrix(x, ncol=length(var.names), byrow=FALSE)[new.idx,]),
                              x.order=x.order,
                              offset=offset[new.idx],
                              Misc=Misc, w=w[new.idx] )
      }
      else
      {
         gbm.obj$data <- list(y=y,x=x,x.order=x.order,offset=offset,Misc=Misc,w=w)
      }
   }
   else
   {
      gbm.obj$data <- NULL
   }

   class(gbm.obj) <- "gbm"
   return(gbm.obj)
}
