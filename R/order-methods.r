# Ordering functions 
# A series of functions for ordering gbm_data and
# training_params objects according to observation id and groupings.
# The latter is only relevant for the Pairwise distribution.

reorder_fit <- function(gbm_fit, distribution_obj) {
  UseMethod("reorder_fit", distribution_obj)
  return(gbm_fit)
}

reorder_fit.CoxPHGBMDist <- function(gbm_fit, distribution_obj) {
  gbm_fit$fit[distribution_obj$time_order] <- gbm_fit$fit
  return(gbm_fit)
}

reorder_fit.PairwiseGBMDist <- function(gbm_fit, distribution_obj) {
  gbm_fit$fit <- gbm_fit$fit[distribution_obj$group_order]
  return(gbm_fit)
}

order_data <- function(gbm_data_obj, distribution_obj, train_params) {
  gbm_data_obj <- order_by_groupings(gbm_data_obj, distribution_obj)
  gbm_data_obj <- order_by_id(gbm_data_obj, train_params)
  gbm_data_obj <- predictor_order(gbm_data_obj, train_params)
  return(gbm_data_obj)
}



order_by_id <- function(gbm_data_obj, train_params) {
  # Check if gbm_data_obj
  check_if_gbm_data(gbm_data_obj)
  gbm_data_obj$x <- gbm_data_obj$x[train_params$id, , drop=FALSE]
  gbm_data_obj$y <- gbm_data_obj$y[train_params$id]
  
  return(gbm_data_obj)
}

order_by_groupings <- function(gbm_data_obj, distribution_obj) {
  # Check if gbm_data_obj
  check_if_gbm_data(gbm_data_obj)
  
  # Check if GBMDist obj
  check_if_gbm_dist(distribution_obj)
  if(!is.null(distribution_obj$group)) {
    gbm_data_obj$y            <- gbm_data_obj$y[distribution_obj$group_order]
    gbm_data_obj$x            <- gbm_data_obj$x[distribution_obj$group_order,,drop=FALSE]
    gbm_data_obj$weights      <- gbm_data_obj$weights[distribution_obj$group_order]
  }
  return(gbm_data_obj)
}