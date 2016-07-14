#' GBM CV errors
#' 
#' Calculate cross validation errors
#' 
#' @usage gbm_cv_errors(gbm_cv_fit, cv_folds, cv_group)
#' 
#' @param gbm_cv_fit a GBMCVFit object containing CV gbm models
#' 
#' 
#' @param cv_folds a positive integer specifying the number of folds to be used in cross-validation of the gbm
#' fit.
#' 
#' @param cv_group vector of integers specifying which row of data belongs to which cv_fold.

gbm_cv_errors <- function(gbm_cv_fit, cv_folds, cv_group) {
  UseMethod("gbm_cv_errors", gbm_cv_fit)
}

gbm_cv_errors.GBMCVFit <- function(gbm_cv_fit, cv_folds, cv_group) {
  in_group <- tabulate(cv_group, nbins=cv_folds)
  cv_error <- vapply(1:cv_folds,
                     function(index) {
                       model <- gbm_cv_fit[[index+1]]
                       model$valid.error * in_group[[index]]
                     }, double(gbm_cv_fit[[1]]$params$num_trees))
  ## this is now a (num_trees, cv_folds) matrix
  ## and now a n.trees vector
  return(rowSums(cv_error) / gbm_cv_fit[[1]]$params$num_train)
}