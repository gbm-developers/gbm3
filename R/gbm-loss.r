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

loss <- function(y, predictions, weights, offset, distribution_obj, baseline=rep(0, length(y))) {
  # Check inputs
  check_if_gbm_dist(distribution_obj)
  if(any(y != as.double(y)) || is.na(y)) {
    stop("Responses must be doubles")
  }
  if(any(predictions != as.double(predictions)) || is.na(predictions)) {
    stop("Predictions must be doubles")
  }
  if(any(weights != as.double(weights)) || is.na(weights)) {
    stop("Weights must be doubles")
  }
  if(any(baseline != as.double(baseline)) || is.na(baseline)) {
    stop("Baseline must consist of doubles")
  }
  
  if((length(y) != length(predictions)) || (length(predictions) != length(wieghts))
      || (length(y) != length(baseline))) {
    stop("Predictions, responses, weightsand baseline all must have the same number
         of elements")
  }
  
  # Add on offset
  if (!is.na(offset) && any(offset != as.double(offset)))
  {
    if((length(offset) != length(predictions))) {
      stop("Offset must be the same length as the prediction ")
    }
    predictions <- offset+prediction
  }
  
  UseMethod("loss", distribution_obj)
}

loss.default <- function(y, predictions, weights, offset, distribution_obj, baseline) {
  stop("loss function not specified for distribution object provided.")
}

loss.AdaBoostGBMDist <- function(y, predictions, weights, offset, distribution_obj, baseline) {
  return(weighted.mean(exp(-(2*y-1)*predictions), weights) - baseline)
}

loss.BernoulliGBMDist <- function(y, predictions, weights, offset, distribution_obj, baseline) {
  return(-2*weighted.mean(y*predictions - log(1+exp(predictions)), weights) - baseline)
}

loss.CoxPHGBMDist <- function(y, predictions, weights, offset, distribution_obj, baseline) {
  stop("Loss method for ",  class(distribution_obj)[1]," not yet supported.")
}

loss.GammaGBMDist <- function(y, predictions, weights, offset, distribution_obj, baseline) {
  stop("Loss method for ",  class(distribution_obj)[1]," not yet supported.")
}

loss.GaussianGBMDist <- function(y, predictions, weights, offset, distribution_obj, baseline) {
  return(weighted.mean((y - predictions)^2, weights) - baseline)
}

loss.HuberizedGBMDist <- function(y, predictions, weights, offset, distribution_obj, baseline) {
  stop("Loss method for ",  class(distribution_obj)[1]," not yet supported.")
}

loss.LaplaceGBMDist <- function(y, predictions, weights, offset, distribution_obj, baseline) {
  return(weighted.mean(abs(y-predictions), weights) - baseline)
}

loss.PairwiseGBMDist <- function(y, predictions, weights, offset, distribution_obj, baseline) {
  return((1 - perf.pairwise(y, predictions, distribution_obj$group, distribution_obj$metric, 
                            weights, distribution_obj$max.rank)) - baseline)
}

loss.PoissonGBMDist <- function(y, predictions, weights, offset, distribution_obj, baseline) {
  return(-2*weighted.mean(y*predictions-exp(predictions), weights) - baseline)
}

loss.QuantileGBMDist <- function(y, predictions, weights, offset, distribution_obj, baseline) {
  stop("Loss method for ",  class(distribution_obj)[1]," not yet supported.")
  
}

loss.TDistGBMDist <- function(y, predictions, weights, offset, distribution_obj, baseline) {
  stop("Loss method for ",  class(distribution_obj)[1]," not yet supported.")
  
}

loss.TweedieGBMDist <- function(y, predictions, weights, offset, distribution_obj, baseline) {
  stop("Loss method for ",  class(distribution_obj)[1]," not yet supported.")
  
}
