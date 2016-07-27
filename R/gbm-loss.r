#' GBMFit loss
#' 
#' S3 method for calculating the loss given data and a GBMDist object
#' 
#' @usage loss(y, predictions, weights, offset, distribution_obj, baseline=rep(0, length(y)))
#' 
#' @param y a vector of responses 
#' 
#' @param predictions vector of predicted respones
#' 
#' @param weights weightings of each point in loss calculation
#' 
#' @param offset offset for each prediction
#' 
#' @param distribution_obj a GBMDist object which determines how the loss will be calculated
#' 
#' @param baseline a vector of doubles specifying the baseline from which the loss is calculated.  This
#' defaults to 0.
#' 
#' @return the loss associated with the fit and distribution - vector of doubles
#' 
#' @export 
loss <- function(y, predictions, weights, offset, distribution_obj, baseline=rep(0, length(y))) {
  # Check inputs
  check_if_gbm_dist(distribution_obj)
  if(is.logical(y) || is.character(y) || any(y != as.double(y)) || is.na(y)) {
    stop("Responses must be doubles")
  }
  if(is.logical(predictions) || is.character(predictions) || any(predictions != as.double(predictions)) || is.na(predictions)) {
    stop("Predictions must be doubles")
  }
  if(is.logical(weights) || is.character(weights) || any(weights != as.double(weights)) || is.na(weights)) {
    stop("Weights must be doubles")
  }
  if(is.logical(baseline) || is.character(baseline) || any(baseline != as.double(baseline)) || is.na(baseline)) {
    stop("Baseline must consist of doubles")
  }
  
  if(is.logical(offset) || is.character(offset) || any(offset != as.double(offset)) || is.na(offset)) {
    stop("Offset must consist of doubles")
  }
  
  if((length(y) != length(predictions)) || (length(predictions) != length(weights))
      || (length(y) != length(baseline))) {
    stop("Predictions, responses, weights and baseline all must have the same number
         of elements")
  }
  
  if((length(offset) != length(predictions))) {
    stop("Offset must be the same length as the prediction ")
  }
   
  UseMethod("loss", distribution_obj)
}

#' @name loss
#' @export
loss.default <- function(y, predictions, weights, offset, distribution_obj, baseline) {
  stop("loss function not specified for distribution object provided.")
}

#' @name loss
#' @export
loss.AdaBoostGBMDist <- function(y, predictions, weights, offset, distribution_obj, baseline) {
  return(weighted.mean(exp(-(2*y-1)*(predictions+offset)), weights) - baseline)
}

#' @name loss
#' @export
loss.BernoulliGBMDist <- function(y, predictions, weights, offset, distribution_obj, baseline) {
  return(-2*weighted.mean(y*(predictions+offset) - log(1+exp(predictions+offset)), weights) - baseline)
}

#' @name loss
#' @export
loss.CoxPHGBMDist <- function(y, predictions, weights, offset, distribution_obj, baseline) {
  stop("Loss method for ",  class(distribution_obj)[1]," not yet supported.")
}

#' @name loss
#' @export
loss.GammaGBMDist <- function(y, predictions, weights, offset, distribution_obj, baseline) {
  stop("Loss method for ",  class(distribution_obj)[1]," not yet supported.")
}

#' @name loss
#' @export
loss.GaussianGBMDist <- function(y, predictions, weights, offset, distribution_obj, baseline) {
  return(weighted.mean((y - predictions - offset)^2, weights) - baseline)
}

#' @name loss
#' @export
loss.HuberizedGBMDist <- function(y, predictions, weights, offset, distribution_obj, baseline) {
  stop("Loss method for ",  class(distribution_obj)[1]," not yet supported.")
}

#' @name loss
#' @export
loss.LaplaceGBMDist <- function(y, predictions, weights, offset, distribution_obj, baseline) {
  return(weighted.mean(abs(y-predictions - offset), weights) - baseline)
}

#' @name loss
#' @export
loss.PairwiseGBMDist <- function(y, predictions, weights, offset, distribution_obj, baseline) {
  return((1 - perf_pairwise(y, predictions+offset, distribution_obj$group, distribution_obj$metric, 
                            weights, distribution_obj$max_rank)) - baseline)
}

#' @name loss
#' @export
loss.PoissonGBMDist <- function(y, predictions, weights, offset, distribution_obj, baseline) {
  return(-2*weighted.mean(y*(predictions+offset)-exp(predictions+offset), weights) - baseline)
}

#' @name loss
#' @export
loss.QuantileGBMDist <- function(y, predictions, weights, offset, distribution_obj, baseline) {
  stop("Loss method for ",  class(distribution_obj)[1]," not yet supported.")
  
}

#' @name loss
#' @export
loss.TDistGBMDist <- function(y, predictions, weights, offset, distribution_obj, baseline) {
  stop("Loss method for ",  class(distribution_obj)[1]," not yet supported.")
  
}

#' @name loss
#' @export
loss.TweedieGBMDist <- function(y, predictions, weights, offset, distribution_obj, baseline) {
  stop("Loss method for ",  class(distribution_obj)[1]," not yet supported.")
  
}
