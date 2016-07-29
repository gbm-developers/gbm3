#' GBM data object
#' Object storing the data used in training a gbm model. 
#' 
#' @usage gbm_data(x, y, weights, offset)
#' 
#' @param x matrix or data-frame of the predictor variables 
#' 
#' @param y matrix or vector of response variables
#' 
#' @param weights vector containing the weights applied to each data row
#' 
#' @param offset vector containing the response variables offsets
#' 
#' @author James Hickey
#' 
#' @return a gbm_data object 
#' 
#' @export gbm_data

gbm_data <- function(x, y, weights, offset) {
    
    # Check inputs 
    verify_data(x, y)
    
    # Convert y if 2-level factor
    y <- convertY(y)
  
    # Check weights and offsets are doubles  
    weights <- check_weights(weights, length(y))
    if(is.null(weights) || is.infinite(weights) || !is.atomic(weights)
       || !is.double(weights)) {
      stop("Weights must be a vector of doubles")
    }
    
    if(is.null(offset) || is.infinite(offset) || !is.atomic(offset) ||
       !is.double(offset) || (length(offset) == 1) ) {
      stop("Offsets must be a vector of doubles")
    }
    
    # Store original data - before any formatting
    data <- cbind(as.data.frame(y), as.data.frame(x))
    
    return(structure(list(x=x, y=y, original_data=data, weights=weights, offset=offset),
                     class="GBMData"))
}