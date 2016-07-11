#' Predictor Order
#' 
#' Function that calculates the predictor variable order for GBM
#' 
#' @usage predictor_order(gbm_data_obj, train_params)
#' 
#' @param gbm_data_obj a \code{GBMData} object.
#' 
#' @param train_params a \code{GBMTrainParams} object.
#'
#' @return an updated gbm_data_object that now contains the predictor variable order.
#' 
#' @export

predictor_order <- function(gbm_data_obj, train_params) {
  gbm_data_obj$x_order <- apply(gbm_data_obj$x[seq_len(train_params$num_train),,drop=FALSE],2,order,na.last=FALSE)-1
  return(gbm_data_obj)
}