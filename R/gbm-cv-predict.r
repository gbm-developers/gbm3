#' Predictions for CV fitted GBM models
#' 
#' Calculate predictions from fitted cv gbm results object
#' 
#' @usage predict(gbm_cv_fit, gbm_data_obj, cv_folds, cv_group, best_iter_cv)
#' 
#' @param gbm_cv_fit a GBMCVFit object containing CV gbm models
#' 
#' @param gbm_data_obj a GBMData object containing all of the data used to fit a gbm model. 
#' 
#' @param cv_folds a positive integer specifying the number of folds to be used in cross-validation of the gbm
#' fit.
#' 
#' @param cv_group vector of integers specifying which row of data belongs to which cv_fold.
#' 
#' @param best_iter_cv number of trees with the smallest cv error

predict.GBMCVFit <- function(gbm_cv_fit, gbm_data_obj, cv_folds, cv_group, best_iter_cv) {
  
  # Extractg data
  data <- gbm_data_obj[seq_len(gbm_cv_fit[[1]]$train_params$num_train), ,drop=FALSE]
  
  ## test cv_group and data match
  if (nrow(data) != length(cv_group)) {
    stop("mismatch between data and cv_group")
  }
  
  num_cols <- 1
  result <- matrix(nrow=nrow(data), ncol=num.cols)
  data_names <- names(data)
  for (ind in seq_len(cv_folds)) {
    flag <- cv_group == ind
    model <- gbm_cv_fit[[ind+1]]
    my_data  <- data[flag, model$variables$var_names, drop=FALSE]
    predictions <- predict(model, newdata=my_data, num_trees=best_iter_cv)
    predictions <- matrix(predictions, ncol=num_cols)
    result[flag,] <- predictions
  }
  
  return(as.numeric(result))
}