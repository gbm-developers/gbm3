#' Convert factors
#' 
#' Function that converts factors or ordered predictors
#' to an appropriate numeric type. 
#' 
#' @usage convert_factor(gbm_data_obj)
#'
#' @param gbm_data_obj a gbm_data object
#'   
#' @return gbm_data object with ordered/factor predictor variables appropriately converted
#' 

convert_factors <- function(gbm_data_obj) {
  check_if_gbm_data(gbm_data_obj)
  gbm_data_obj$x <- sapply(gbm_data_obj$x, function(col) if(is.ordered(col) || is.factor(col)){
                                                         col <- as.numeric(factor(col))-1 } )
  return(gbm_data_obj)
}