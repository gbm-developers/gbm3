#' @describeIn gbm add new trees to a gbm object
#' @export
gbm.more <- function(object,
                     n.new.trees = 100,
                     data = NULL,
                     weights = NULL,
                     offset = NULL,
                     verbose = NULL,
                     n.threads = 1)
{
   theCall <- match.call()
   
   # Get 
   nTrainRows  <- object$nTrain
   nTrainPats <- object$nTrainPats
   patient.id <- object$patient.id
   StrataVec <- as.vector(object$strata)
   sortedVec <- as.vector(object$sorted)
   prior.node.coeff.var <- as.double(object$prior.node.coeff.var)
   
   if (object$distribution$name != "pairwise") {
      distribution.call.name <- object$distribution$name
   } else {
     distribution.call.name <- sprintf("pairwise_%s", object$distribution$metric)
   }

   if(is.null(object$Terms) && is.null(object$data)) {
     stop("The gbm model was fit using gbm.fit (rather than gbm) and keep.data was set to FALSE. gbm.more cannot locate the dataset.")
   } else if(is.null(object$data) && is.null(data)) {
     stop("keep.data was set to FALSE on original gbm call and argument 'data' is NULL")
   }  else if(is.null(object$data)) {
     
      m <- eval(object$m, parent.frame())

      Terms <- attr(m, "terms")
      a <- attributes(Terms)

      y <- as.vector(model.extract(m, "response"))
      offset <- model.extract(m,offset)
      x <- model.frame(delete.response(Terms),
                       data,
                       na.action=na.pass)

      w <- weights
      if(length(w)==0) w <- rep(1, nrow(x))

      if (object$distribution$name != "pairwise")
      {
         w <- w*length(w)/sum(w) # normalize to N
      }

      if(is.null(offset) || (offset==0))
      {
         offset <- NA
      }
      Misc <- NA

      if(object$distribution$name == "coxph") 
      {
        if(class(y)!="Surv")
        {
          stop("Outcome must be a survival object Surv(time1, failure) or Surv(time1, time2, failure)")
        }
        
        
        # Patients are split into train and test, and are ordered by
        # strata
        # Define number of tests
        n.test <- cRows - nTrainRows
        y <- as.numeric(y)[1:cRows]
        
        # Sort response according to strata
        # i.order sets order of outputs
        if (attr(y, "type") == "right")
        {
          sorted <- c(order(StrataVec[1:nTrainRows], -y[1:nTrainRows, 1]), order(StrataVec[(nTrainRows+1):cRows], -y[(nTrainRows+1):cRows, 1])) 
          i.order <- c(order(StrataVec[1:nTrainRows], -y[1:nTrainRows, 1]), order(StrataVec[(nTrainRows+1):cRows], -y[(nTrainRows+1):cRows, 1]) + nTrainRows)
        }
        else if (attr(y, "type") == "counting") 
        {
          sorted <- cbind(c(order(StrataVec[1:nTrainRows], -y[1:nTrainRows, 1]), order(StrataVec[(nTrainRows+1):cRows], -y[(nTrainRows+1):cRows, 1])),
                          c(order(StrataVec[1:nTrainRows], -y[1:nTrainRows, 2]), order(StrataVec[(nTrainRows+1):cRows], -y[(nTrainRows+1):cRows, 2])))
          i.order <- c(order(StrataVec[1:nTrainRows], -y[1:nTrainRows, 1]), order(StrataVec[(nTrainRows+1):cRows], -y[(nTrainRows+1):cRows, 1]) + nTrainRows)
        }
        else
        {
          stop("Survival object must be either right or counting type.")
        }
        
        # Add in sorted column and strata
        sortedVec <- sorted-1L
        
        # Set ties here for the moment
        Misc <- list("ties" = "efron")
        
      } else if(object$distribution$name == "tdist" ){
        Misc <- object$distribution$df
      } else if (object$distribution$name == "pairwise"){

         # Check if group names are valid
         distribution.group <- object$distribution$group
         i <- match(distribution.group, colnames(data))
         if (any(is.na(i))) {
            stop("Group column does not occur in data: ", distribution.group[is.na(i)])
         }

         # construct group index
         group <- factor(do.call(paste, c(data[,distribution.group, drop=FALSE], sep=":")))

         # Check that weights are constant across groups
         if ((!missing(weights)) && (!is.null(weights))) {
            w.min <- tapply(w, INDEX=group, FUN=min)
            w.max <- tapply(w, INDEX=group, FUN=max)

            if (any(w.min != w.max)) {
               stop("For distribution 'pairwise', all instances for the same group must have the same weight")
            }

            # Normalize across groups
            w <- w * length(w.min) / sum(w.min)
         }

         # Shuffle groups, to remove bias when splitting into train/test set and/or CV folds
         perm.levels  <- levels(group)[sample(1:nlevels(group))]
         group        <- factor(group, levels=perm.levels)

         # The C function expects instances to be sorted by group and descending by target

         ord.group    <- object$ord.group
         group        <- group[ord.group]
         y            <- y[ord.group]
         x            <- x[ord.group,,drop=FALSE]
         w            <- x[ord.group]
         object$fit   <- object$fit[ord.group] # object$fit is stored in the original order

         # Split into train and validation set, at group boundary
         num.groups.train <- max(1, round(object$train.fraction * nlevels(group)))

         # include all groups up to the num.groups.train
         nTrainRows           <- max(which(group==levels(group)[num.groups.train]))

         metric <- object$distribution[["metric"]]

         if (is.element(metric, c("mrr", "map")) && (!all(is.element(y, 0:1)))) {
            stop("Metrics 'map' and 'mrr' require the response to be in {0,1}")
         }

         # Cut-off rank for metrics
         # We pass this argument as the last element in the Misc vector
         # Default of 0 means no cutoff

         max.rank <- 0
         if (!is.null(object$distribution[["max.rank"]]) && object$distribution[["max.rank"]] > 0) {
            if (is.element(metric, c("ndcg", "mrr"))) {
               max.rank <- object$distribution[["max.rank"]]
            } else {
               stop("Parameter 'max.rank' cannot be specified for metric '", metric, "', only supported for 'ndcg' and 'mrr'")
            }
         }

         Misc <- c(group, max.rank)

      }

      # create index upfront... subtract one for 0 based order
      x.order <- apply(x[1:nTrainRows,,drop=FALSE],2,order,na.last=FALSE)-1
      x <- data.matrix(x)
      cRows <- nrow(x)
      cCols <- ncol(x)
   } else {
      y       <- object$data$y
      x       <- object$data$x
      x.order <- object$data$x.order
      offset  <- object$data$offset
      Misc    <- object$data$Misc
      w       <- object$data$w
      nTrainRows  <- object$nTrain
      nTrainPats <- object$nTrainPats
      cRows   <- length(y)
      cCols   <- length(x)/cRows
      
     
      if (object$distribution$name == "pairwise") 
      {
         object$fit   <- object$fit[object$ord.group] # object$fit is stored in the original order
      }
   }

   # Get size of Response frame
   y <- as.vector(y)
   cRowsY <- nrow(y)
   cColsY <- ncol(y)
   
   if(is.null(cRowsY))
   {
     cRowsY <- length(y)
     cColsY <- 1
   }
   
   # Make sorted vec into a matrix
   if(cColsY > 2)
   {
     cRowsSort <- dim(sortedVec)[1]
     cColsSort <- dim(sortedVec)[2]
   }
   else
   {
     cRowsSort <- length(sortedVec)
     cColsSort <- 1
   }
   

   
   if(is.null(verbose)) {
      verbose <- object$verbose
   }
   x <- as.vector(x)

   gbm.obj <- .Call("gbm",
                    Y = matrix(y, cRowsY, cColsY),
                    Offset = as.double(offset),
                    X = matrix(x, cRows, cCols),
                    X.order = as.integer(x.order),
                    sorted=matrix(sortedVec, cRowsSort, cColsSort),
                    Strata = as.integer(StrataVec),
                    weights = as.double(w),
                    Misc = as.double(Misc),
                    prior.node.coeff.var = as.double(prior.node.coeff.var),
                    patient.id = as.integer(patient.id),
                    var.type = as.integer(object$var.type),
                    var.monotone = as.integer(object$var.monotone),
                    distribution = as.character(distribution.call.name),
                    n.trees = as.integer(n.new.trees),
                    interaction.depth = as.integer(object$interaction.depth),
                    n.minobsinnode = as.integer(object$n.minobsinnode),
                    shrinkage = as.double(object$shrinkage),
                    bag.fraction = as.double(object$bag.fraction),
                    nTrain = as.integer(nTrainRows), #Should this be as.double(train.fraction) ? 
                    nTrainPats = as.integer(nTrainPats),
                    mFeatures = as.integer(object$mFeatures),
                    fit.old = as.double(object$fit),
                    n.cat.splits.old = as.integer(length(object$c.splits)),
                    n.trees.old = as.integer(object$n.trees),
                    n.threads,
                    verbose = as.integer(verbose),
                    PACKAGE = "gbm")
   
   gbm.obj$initF         <- object$initF
   gbm.obj$train.error   <- c(object$train.error, gbm.obj$train.error)
   gbm.obj$valid.error   <- c(object$valid.error, gbm.obj$valid.error)
   gbm.obj$oobag.improve <- c(object$oobag.improve, gbm.obj$oobag.improve)
   gbm.obj$trees         <- c(object$trees, gbm.obj$trees)
   gbm.obj$c.splits      <- c(object$c.splits, gbm.obj$c.splits)

   # cv.error not updated when using gbm.more
   gbm.obj$cv.error      <- object$cv.error
   gbm.obj$cv.folds      <- object$cv.folds

   gbm.obj$n.trees        <- length(gbm.obj$trees)
   gbm.obj$distribution   <- object$distribution
   gbm.obj$train.fraction <- object$train.fraction
   gbm.obj$shrinkage      <- object$shrinkage
   gbm.obj$bag.fraction   <- object$bag.fraction
   gbm.obj$var.type       <- object$var.type
   gbm.obj$var.monotone   <- object$var.monotone
   gbm.obj$var.names      <- object$var.names
   gbm.obj$interaction.depth <- object$interaction.depth
   gbm.obj$n.minobsinnode    <- object$n.minobsinnode
   gbm.obj$num.classes       <- object$num.classes
   gbm.obj$nTrain            <- object$nTrain
   gbm.obj$nTrainPats <- object$nTrainPats
   gbm.obj$mFeatures         <- object$mFeatures
   gbm.obj$response.name     <- object$response.name
   gbm.obj$Terms             <- object$Terms
   gbm.obj$var.levels        <- object$var.levels
   gbm.obj$verbose           <- verbose
   gbm.obj$prior.node.coeff.var <- object$prior.node.coeff.var
   gbm.obj$patient.id <- object$patient.id
   gbm.obj$strata <- object$strata
   gbm.obj$sorted <- object$sorted

   if(object$distribution$name == "coxph")
   {
      gbm.obj$fit[i.order] <- gbm.obj$fit
   }

   if (object$distribution$name == "pairwise") {
      # Data has been reordered according to queries.
      # We need to permute the fitted values to correspond
      # to the original order.
      gbm.obj$fit <- gbm.obj$fit[order(object$ord.group)]
      object$fit  <- object$fit[order(object$ord.group)]
      gbm.obj$ord.group <- object$ord.group
   }

   if(!is.null(object$data)) {
      gbm.obj$data <- object$data
   } else {
      gbm.obj$data <- NULL
   }
   gbm.obj$m <- object$m
   gbm.obj$call <- theCall

   class(gbm.obj) <- "gbm"
   return(gbm.obj)
}
