#' GBM data object
#' Object storing the data used in training a gbm model. 
#' 
#' @usage gbm_data(x, y, x_order, weights, offset, id_order)
#' 
#' @param x matrix or data-frame of the predictor variables 
#' 
#' @param y matrix or vector of response variables
#' 
#' @param weights vector containing the weights applied to each data row
#' 
#' @param offset vector containing the response variables offsets
#' 
#' @param id_order the ordering of the variables according to the observation id.
#' 
#' 
#' @return a gbm_data object 
#' 
#' @export gbm_data

gbm_data <- function(x, y, weights, offset, id_order, distribution_obj) {
    
    # Check inputs 
    verify_data(x, y)
    
    # Check weights and offsets are doubles  
    if(is.null(weights) && is.infinite(weights) && !is.atomic(weights)
       && !is.double(weights)) {
      stop("Weights must be a vector of doubles")
    }
    
    if(is.null(offset) && is.infinite(offset) && !is.atomic(offset)
       && !is.double(offset)) {
      stop("Offsets must be a vector of doubles")
    }
    
    # Normalise Weights
    if(!(distribution_obj$name %in% available_distributions())) {
      stop("Distribution not recognised - see available_distributions for
           supported distributions")
    }
    w <- checkWeights(w, length(y))
    if(distribution_obj$name != "Pairwise") {
      w <- w*length(w)/sum(w)
      message("Distribution selected is not Pairwise - normalizing weights")
    } 
  

    offset <- checkOffset(offset, y, distribution)
    # Check id_order
    if(!is.atomic(id_order) || any(id_order != trunc(id_order)) ||is.infinite(id_order)
       || is.null(id_order) || (length(id_order) != length(y)) || (any(id_order) < 0)) {
      stop("Ordering of observations must be a vector of positive whole numbers
           of the same length as the number of rows of data.")
    }
    
    return(structure(list(x=x, y=y, weights=weights, offset=offset, id_order=id_order),
                     class="GBMData"))
}