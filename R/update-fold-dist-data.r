#' Recreate distribution data for CV folds
#' 
#' Modify distribution object data for CV gbm fit.
#' 
#' @usage update_fold_dist_data(gbm_dist_obj, gbm_data_obj, train_params, rows_in_training, rows_in_fold)
#'  
#' @param gbm_dist_obj a GBMDist object properly formatted for gbm_fit.
#'
#' @param gbm_data_obj  a GBMData object containing correctly ordered and validated data.
#' 
#' @param train_params a validated GBMTrainParams object
#' 
#' @param rows_in_training vector of logicals that determine what data rows are in the training set.
#' 
#' @param rows_in_fold vector of logicals indicating whether a row of training data is in the fold or not.
#' 
#' @return modified \code{GBMDist} object for CV fit - strata and groupings updated if needed.
#' 
#' @export update_fold_dist_data update_fold_dist_data.default update_fold_dist_data.CoxPHGBMDist update_fold_dist_data.PairwiseGBMDist

update_fold_dist_data <- function(gbm_dist_obj, gbm_data_obj, train_params, rows_in_training, rows_in_fold) {
  check_if_gbm_dist(gbm_dist_obj)
  check_if_gbm_data(gbm_data_obj)
  check_if_gbm_train_params(train_params)
  if(!is.atomic(rows_in_fold) || any(!is.logical(rows_in_fold)) ||
     (length(rows_in_fold[rows_in_fold == FALSE]) != train_params$num_train)) {
    stop("rows_in_fold must be a vector of logicals of length the number of training rows")
  }
  UseMethod("update_fold_dist_data", gbm_dist_obj)
}

update_fold_dist_data.default <- function(gbm_dist_obj, gbm_data_obj, train_params, rows_in_training, rows_in_fold) {
  return(gbm_dist_obj)
}

update_fold_dist_data.CoxPHGBMDist <- function(gbm_dist_obj, gbm_data_obj, train_params, rows_in_training, rows_in_fold) {
  # Reset strata using folds
  gbm_dist_obj$strata <- NULL
  gbm_dist_obj$sorted <- NULL
  gbm_dist_obj$time_order <- NULL
  gbm_dist_obj <- create_strata(gbm_data_obj, train_params, gbm_dist_obj)
  return(gbm_dist_obj)
}

update_fold_dist_data.PairwiseGBMDist <- function(gbm_dist_obj, gbm_data_obj, train_params, rows_in_training, rows_in_fold) {
  gbm_dist_obj$group <- c(gbm_dist_obj$group[rows_in_training][!rows_in_fold],
                          gbm_dist_obj$group[rows_in_training][rows_in_fold])
  return(gbm_dist_obj)
}
