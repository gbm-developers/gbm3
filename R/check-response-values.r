# Check Response Values
#
# check responses are consistent with the distribution specified.
# 
# @usage check_response_values(y, distribution_obj)
# 
# @param y a matrix or vector of outcomes
# 
# @param distribution_obj a GBMDist object
# 
# @author James Hickey
#
# @return Warnings/Errors
# 

check_response_values <- function(distribution_obj, y) {
  # Check inputs
  if( !(class(distribution_obj)[1] %in%
        paste0(available_distributions(),"GBMDist")) ) {
    stop("Distribution object not recognised - use gbm_dist to create a valid object")
  }
  
  if(!is.matrix(y) && !is.data.frame(y) &&
     !is.atomic(y)) {
    stop("Responses must be in a dataframe, matrix or vector")
  } 
  
  # Call correct method
  UseMethod("check_response_values", distribution_obj)
  
}

check_response_values.AdaBoostGBMDist <-function(distribution_obj, y) {
  if(!all(is.element(y,0:1))) {
    stop("This version of AdaBoost requires the response to be in {0,1}")
  }
}

check_response_values.BernoulliGBMDist <-function(distribution_obj, y) {
  if(!all(is.element(y, 0:1))) {
    stop("Bernoulli requires the response to be in {0,1}")
  }
}

check_response_values.CoxPHGBMDist <-function(distribution_obj, y) {
  if(!inherits(y, "Surv")) {
    stop("Outcome must be a survival object Surv(time1, failure) 
         or Surv(time1, time2, failure)")
  }
  
  # Check length if not default
  if(!is.na(distribution_obj$original_strata_id) && 
     (length(distribution_obj$original_strata_id) != nrow(y)) ){
    stop("Strata indices must be provided for every data point")
  }
}

check_response_values.GammaGBMDist <-function(distribution_obj, y) {
  if(any(y<0)) {
    stop("Gamma requires the response to be positive")
  }
}

check_response_values.GaussianGBMDist <-function(distribution_obj, y) {}

check_response_values.HuberizedGBMDist <-function(distribution_obj, y) {
  if(!all(is.element(y,0:1))) {
    stop("Huberized square hinged loss requires the response to be in {0,1}")
  }
}

check_response_values.LaplaceGBMDist <-function(distribution_obj, y) {}

check_response_values.PairwiseGBMDist <-function(distribution_obj, y) {
  if (any(y<0)) {
    stop("targets for 'pairwise' should be non-negative")
  }
  
  if (is.element(distribution_obj$metric, c("mrr", "map")) && (!all(is.element(y, 0:1)))) {
    stop("Metrics 'map' and 'mrr' require the response to be in {0,1}")
  }
}

check_response_values.PoissonGBMDist <-function(distribution_obj, y) {
  if(any(y != trunc(y))) {
    stop("Poisson requires the response to be a positive integer")
  }
}

check_response_values.QuantileGBMDist <-function(distribution_obj, y) {}

check_response_values.TDistGBMDist <-function(distribution_obj, y) {}

check_response_values.TweedieGBMDist <-function(distribution_obj, y) {
  if(any(y<0)) {
    stop("Tweedie requires the response to be positive")
  }
}
