# Cross-validate a gbm
# 
# Functions for cross-validating gbm. These functions are used internally and
# are not intended for end-user direct usage.
# 
# These functions are not intended for end-user direct usage, but are used
# internally by \code{gbm}.
# 
# aliases: gbmCrossVal gbmCrossValModelBuild gbmDoFold gbmCrossValErr
# gbmCrossValPredictions
gbmCrossVal <- function(cv.folds, nTrain, n.cores,
                        class.stratify.cv, data,
                        x, y, offset, distribution, w, var.monotone,
                        n.trees, interaction.depth, n.minobsinnode,
                        shrinkage, bag.fraction, mFeatures,
                        var.names, response.name, group, lVerbose, keep.data,
                        fold.id, tied.times.method) {
  i.train <- 1:nTrain
  cv.group <- getCVgroup(distribution$name, class.stratify.cv, y,
                         i.train, cv.folds, group, fold.id)
  ## build the models
  cv.models <- gbmCrossValModelBuild(cv.folds, cv.group, n.cores,
                                     i.train, x, y, offset,
                                     distribution, w, var.monotone,
                                     n.trees, interaction.depth,
                                     n.minobsinnode, shrinkage,
                                     bag.fraction, mFeatures, var.names,
                                     response.name, group, lVerbose, keep.data, 
                                     nTrain, tied.times.method)

  # First element is final model
  all.model <- cv.models[[1]]
  cv.models <- cv.models[-1]
  ## get the errors
  cv.error  <- gbmCrossValErr(cv.models, cv.folds, cv.group, nTrain, n.trees)
  best.iter.cv <- which.min(cv.error)

  ## get the predictions
  predictions <- gbmCrossValPredictions(cv.models, cv.folds, cv.group,
                                        best.iter.cv, distribution,
                                        data[i.train,,drop=FALSE], y)

  all.model$fold.id <- fold.id

  list(error=cv.error,
       predictions=predictions,
       all.model=all.model)
}

## Get the gbm cross-validation error
gbmCrossValErr <- function(cv.models, cv.folds, cv.group, nTrain, n.trees) {
  in.group <- tabulate(cv.group, nbins=cv.folds)
  cv.error <- vapply(1:cv.folds,
                     function(index) {
                       model <- cv.models[[index]]
                       model$valid.error * in.group[[index]]
                     }, double(n.trees))
  ## this is now a (n.trees, cv.folds) matrix

  ## and now a n.trees vector
  rowSums(cv.error) / nTrain
}

## Get the predictions for GBM cross validation
##
## This function is not as nice as it could be (leakage of y)
gbmCrossValPredictions <- function(cv.models, cv.folds, cv.group,
                                   best.iter.cv, distribution, data, y) {
  ## test cv.group and data match
  if (nrow(data) != length(cv.group)) {
    stop("mismatch between data and cv.group")
  }

  num.cols <- 1
  result <- matrix(nrow=nrow(data), ncol=num.cols)
  ## there's no real reason to do this as other than a for loop
  data.names <- names(data)
  for (ind in 1:cv.folds) {
    ## these are the particular elements
    flag <- cv.group == ind
    model <- cv.models[[ind]]
    ## the %in% here is to handle coxph
    my.data  <- data[flag, model$var.names, drop=FALSE]
    predictions <- predict(model, newdata=my.data, n.trees=best.iter.cv)
    predictions <- matrix(predictions, ncol=num.cols)
    result[flag,] <- predictions
  }

  as.numeric(result)
}


## Perform gbm cross-validation
##
## This function has far too many arguments.
gbmCrossValModelBuild <- function(cv.folds, cv.group, n.cores, i.train,
                                  x, y, offset, distribution,
                                  w, var.monotone, n.trees,
                                  interaction.depth, n.minobsinnode,
                                  shrinkage, bag.fraction, mFeatures,
                                  var.names, response.name,
                                  group, lVerbose, keep.data, nTrain, tied.times.method) {
  ## set up the cluster and add a finalizer
  cluster <- gbmCluster(n.cores)
  on.exit(if (!is.null(cluster)){ parallel::stopCluster(cluster) })

  ## get ourselves some random seeds
  seeds <- as.integer(runif(cv.folds, -(2^31 - 1), 2^31))

  ## now do the cross-validation model builds
  if ( ! is.null(cluster) ){
    parallel::parLapply(cl=cluster, X=0:cv.folds,
            gbmDoFold, i.train, x, y, offset, distribution,
            w, var.monotone, n.trees,
            interaction.depth, n.minobsinnode, shrinkage,
            bag.fraction, mFeatures,
            cv.group, var.names, response.name, group, seeds, lVerbose, keep.data, nTrain, tied.times.method)
  }
  else {
    lapply(X=0:cv.folds,
            gbmDoFold, i.train, x, y, offset, distribution,
            w, var.monotone, n.trees,
            interaction.depth, n.minobsinnode, shrinkage,
            bag.fraction, mFeatures,
            cv.group, var.names, response.name, group, seeds, lVerbose, keep.data, nTrain, tied.times.method)
  }
}
