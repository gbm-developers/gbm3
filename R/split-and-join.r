#' Rearrange training set and distribution data for CV
#' 
#' Methods converts the gbm data and distribution parameters for CV fitting
#' 
#' @usage split_and_join(gbm_data_obj, train_params, rows_in_training, rows_in_fold)
#' 
#' @param gbm_data_obj  a GBMData  containing correctly ordered and validated data.
#' 
#' @param train_params a validated GBMTrainParams object
#' 
#' @param rows_in_training vector of logicals that determine what data rows are in the training set.
#' 
#' @param rows_in_fold vector of logicals indicating whether a row of training data is in the fold or not.
#' 
#' @return gbm_data_obj with validation data moved to end of fields and recalculation of predictor ordering
#' 
#' @export

split_and_join <- function(gbm_data_obj, train_params, rows_in_training, rows_in_fold) {
  require("survival")
  
  # Initial checks
  check_if_gbm_data(gbm_data_obj)
  check_if_gbm_train_params(train_params)
  if(!is.atomic(rows_in_training) || any(!is.logical(rows_in_training))) {
    stop("rows_in_training must be a vector of logicals")
  }
  if(!is.atomic(rows_in_fold) || any(!is.logical(rows_in_fold)) ||
     (length(rows_in_fold[rows_in_fold==FALSE]) != train_params$num_train_rows)) {
    stop("rows_in_fold must be a vector of logicals of length the number of training rows")
  }
  
  # Get Validation and training folds
  gbm_data_obj_train <- gbm_data_obj
  
  # Predictors
  x_valid <- as.data.frame(gbm_data_obj$x[rows_in_training, ,drop=FALSE][rows_in_fold, ,drop=FALSE])
  gbm_data_obj_train$x <- as.data.frame(gbm_data_obj$x[rows_in_training, ,drop=FALSE][!rows_in_fold, ,drop=FALSE])
  
  # Responses 
  y_valid <- as.data.frame(as.matrix(gbm_data_obj$y)[rows_in_training, ,drop=FALSE][rows_in_fold, ,drop=FALSE])
  gbm_data_obj_train$y <- as.data.frame(as.matrix(gbm_data_obj$y)[rows_in_training, ,drop=FALSE][!rows_in_fold, ,drop=FALSE])
  
  # Offset, weights and predictor order for validation
  offset_valid <- gbm_data_obj$offset[rows_in_training][rows_in_fold]
  weights_valid <- gbm_data_obj$weights[rows_in_training][rows_in_fold]
  x_order_valid <- as.data.frame(subset(gbm_data_obj$x_order, rows_in_fold, drop=FALSE))
  gbm_data_obj_train$offset <- gbm_data_obj$offset[rows_in_training][!rows_in_fold]
  gbm_data_obj_train$weights <- gbm_data_obj$weights[rows_in_training][!rows_in_fold]

  # Reorder predictors for fitting
  gbm_data_obj_train <- predictor_order(gbm_data_obj_train, train_params)
  
  # Copy column names - required for recombination
  colnames(y_valid) <- colnames(gbm_data_obj_train$y)
  
  # Recombine - and convert to correct format
  gbm_data_obj$x <- rbind(gbm_data_obj_train$x, x_valid)
  gbm_data_obj$y <- rbind(gbm_data_obj_train$y, y_valid)
  gbm_data_obj$offset <- c(gbm_data_obj_train$offset, offset_valid)
  gbm_data_obj$weights <- c(gbm_data_obj_train$weights, weights_valid)
  gbm_data_obj$x_order <- as.matrix(gbm_data_obj_train$x_order)
  
  return(gbm_data_obj)
}



