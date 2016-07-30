# Ordering methods
# 
# Series of functions that are used to re-order the fitted gbm model,
# depending on what distribution was used in the initial call to \code{\link{gbmt_fit}}.
# Currently reordering only takes place for generalized boosted models produced
# using the CoxPH or Pairwise distributions.  These methods are called internally
# on completion of a call to \code{\link{gbmt_fit}}.
# 
# @usage reorder_fit(gbm_fit_obj_obj, distribution_obj)
# 
# @param gbm_fit_obj_obj a \code{GBMFit} object produced by a previous call to \code{gbmt} or \code{gbmt_fit}.
# 
# @param distributin_obj the \code{GBMDist} object used to produce the gbm_fit_obj_obj.
# 
# @return a \code{GBMFit} object with an appropriately order fit.
#  
# @author James Hickey
#

#### Fit reordering methods ####
reorder_fit <- function(gbm_fit_obj, distribution_obj) {
  check_if_gbm_fit(gbm_fit_obj)
  UseMethod("reorder_fit", distribution_obj)
}

reorder_fit.default <- function(gbm_fit_obj, distribution_obj) {
  return(gbm_fit_obj)
}


reorder_fit.CoxPHGBMDist <- function(gbm_fit_obj, distribution_obj) {
  gbm_fit_obj$fit[distribution_obj$time_order] <- gbm_fit_obj$fit
  return(gbm_fit_obj)
}

reorder_fit.PairwiseGBMDist <- function(gbm_fit_obj, distribution_obj) {
  gbm_fit_obj$fit <- gbm_fit_obj$fit[order(distribution_obj$group_order)]
  return(gbm_fit_obj)
}

#### Data ordering - internal methods ####
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
  if(!is.null(distribution_obj$group_order)) {
    gbm_data_obj$y            <- gbm_data_obj$y[distribution_obj$group_order]
    gbm_data_obj$x            <- gbm_data_obj$x[distribution_obj$group_order,,drop=FALSE]
    gbm_data_obj$weights      <- gbm_data_obj$weights[distribution_obj$group_order]
  }
  return(gbm_data_obj)
}