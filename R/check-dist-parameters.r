# Check Distribution Parameters
# 
# Check if function parameters for creating distribution
# object are of correct form - internal function.
# 
# @usage check_dist_params(empty_obj, ...)
# 
# @param empty_obj A skeleton distribution object 
# (see \code{\link{empty_distribution()}}).
# 
# @param ... Parameters defining how to create a full distribution object,
#  see \code{\link{gbm_dist}}.
#  
# @author James Hickey
#
# @return warnings/Errors
# 


check_dist_params <- function(empty_obj, ...) {
  UseMethod("check_dist_params", empty_obj)
}

check_dist_params.AdaBoostGBMDist <- function(empty_obj, ...) {
  if(length(list(...)) > 0) {
    warning("The ", class(empty_obj)[1], "class does not use any additional
            parameters in construction.")
  }
}

check_dist_params.BernoulliGBMDist <- function(empty_obj, ...) {
  if(length(list(...)) > 0) {
    warning("The ", class(empty_obj)[1], "class does not use any additional
            parameters in construction.")
  }
}

check_dist_params.CoxPHGBMDist <- function(empty_obj, strata, sorted, ties, prior_node_coeff, ...) {
  # Check if additional parameters specified
  if((length(list(...)) > 0)) {
    warning("The ", class(empty_obj)[1], "class only requires 4 additional
            parameters for construction - others provided are ignored.")
  }
  # Check ties
  if(!exists("ties")) {
    stop("Ties parameter not specified - distribution can't be constructed")
  } else if(!is.character(ties)) {
    stop("Ties parameter must be a character vector - distribution
         can't be constructed")
  } else if(!(ties %in% c("breslow", "efron"))) {
    stop("Ties parameter must be either 'breslow' or 'efron'")
  }
  
  # Check strata
  if(!exists("strata")) {
    stop("Strata not specified - distribution could not be constructed")
  } else if (is.null(strata)) {
    stop("Strata should not be NULL")
  } else if(!is.na(strata) && ((!(is.atomic(strata)) || is.infinite(strata)
            || (any(strata != as.integer(strata)) && !all(is.factor(strata))) )) ) {
    stop("Strata parameter must be an atomic of integers or factors")
  } 
  
  # Check sorted
  if(!exists("sorted")) {
    stop("Sorted not specified - distribution could not be constructed")
  } else if(is.null(sorted)) {
    stop("Sorted parameter cannot be NULL")
  } else if(!is.na(sorted) && ((!(is.atomic(sorted)) || is.finite(sorted)
          || !isTRUE(all(sorted == as.integer(sorted)))) ) ) {
    stop("Sorted parameter must be an atomic of integers")
  } 
  
  # Check coeff
  if(!exists("prior_node_coeff")) {
    stop("Prior node coefficient of variation not specified - distribution could
         not be constructed")
  } else if(!is.double(prior_node_coeff) || is.infinite(prior_node_coeff) ||
            (length(prior_node_coeff) > 1)) {
    stop("Prior node coefficient not a finite double - distribution could not be
         constructed")
  } 
}

check_dist_params.GammaGBMDist <- function(empty_obj, ...) {
  if(length(list(...)) > 0) {
    warning("The ", class(empty_obj)[1], "class does not use any additional
            parameters in construction.")
  }
}

check_dist_params.GaussianGBMDist <- function(empty_obj, ...) {
  if(length(list(...)) > 0) {
    warning("The ", class(empty_obj)[1], " class does not use any additional
            parameters in construction.")
  }
}

check_dist_params.HuberizedGBMDist <- function(empty_obj, ...) {
  if(length(list(...)) > 0) {
    warning("The ", class(empty_obj)[1], "class does not use any additional
            parameters in construction.")
  }
}

check_dist_params.LaplaceGBMDist <- function(empty_obj, ...) {
  if(length(list(...)) > 0) {
    warning("The ", class(empty_obj)[1], "class does not use any additional
            parameters in construction.")
  }
}



check_dist_params.PairwiseGBMDist <- function(empty_obj, group, metric,
                                              max_rank, group_index, ...) {
  # Check if parameters are specified
  if(length(list(...)) > 0) {
    warning("The ", class(empty_obj)[1], "class only requires 3 additional
            parameters for construction - others provided are ignored")
  }
  
  # Check group is specified correctly
  if(!exists("group")) {
    stop("Group not specified - ", class(empty_obj)[1], " cannot be initialized.")
  } else if(!is.atomic(group) || (length(group) > 1) || !is.character(group)) {
    stop("Group parameter isn't a character atomic - ", class(empty_obj)[1], 
         " cannot be constructed.")
  }
  
  # Check Metric
  if(!exists("metric")) {
    stop("Metric not specified - ", class(empty_obj)[1], " cannot be constructed.")
  } else if(!(metric %in% c("conc", "mrr", "map", "ndcg"))) {
    stop("Metric specified not supported by gbm package")
  }
    
  # Check max_rank is specified correctly
  if(!exists("max_rank")) {
    stop("Max rank not specified - check default settings")
    
  } else if(!is.double(max_rank) || is.infinite(max_rank) || (length(max_rank) > 1)
            || (max_rank < 0.0)) {
    stop("Max rank provided is not a finite double greater than zero - distribution cannot be constructed")
  } else if((max_rank != 0.0) && (metric %in% c("conc", "map"))){
    stop("Max rank cannot be specified for metrics - conc and map")
  }
  
  # Check group_index
  if(!is.null(group_index)) {
    message("group_index provided this will be updated according to group and data if calling gbmt.")
    if(!is.vector(group_index) || any(group_index != as.integer(group_index))) {
    stop("if group indices are provided it must be a vector of integers")
    }
  }
}

check_dist_params.PoissonGBMDist <- function(empty_obj, ...) {
  if(length(list(...)) > 0) {
    message("The ", class(empty_obj)[1], "class does not use any additional
            parameters in construction.")
  }
}

check_dist_params.QuantileGBMDist <- function(empty_obj, alpha, ...) {
  # Check if parameters are specified
  if(length(list(...)) > 0) {
    warning("The ", class(empty_obj)[1], "class only requires one additional
            parameter for construction - others are ignored")
  }
  
  if(!exists("alpha")) {
    stop("Alpha is not specified - distribution cannot be specified")
  } else if (!is.double(alpha) || (length(alpha) > 1) || is.infinite((alpha))
             || (alpha < 0.0) || (alpha > 1.0)) {
    stop("Alpha provided is not a finite double between 0.0 and 1.0 - distribution 
         cannot be constructed")
  }
}

check_dist_params.TDistGBMDist <- function(empty_obj, df, ...) {
  # Check if parameters are specified
  if(length(list(...)) > 0) {
    warning("The ", class(empty_obj)[1], "class only requires one additional
            parameter for construction - others are ignored")
  }

  if(!exists("df")) {
    stop("Degrees of freedom (df) is not specified - distribution cannot be specified")
  } else if (!is.double(df) || (length(df) > 1) 
             || is.infinite((df)) || df < 0.0) {
    stop("df provided is not a finite double bigger than 0.0 - distribution 
         cannot be constructed")
  }
  
}

check_dist_params.TweedieGBMDist <- function(empty_obj, power, ...) {
  # Check if parameters are specified
  if(length(list(...)) > 0) {
    warning("The ", class(empty_obj)[1], "class only requires one additional
            parameters for construction - others provided are ignored")
  }
  
  if(!exists("power")) {
    stop("Power of distribution (power) is not specified - distribution cannot be specified")
  } else if (!(is.double(power)) || (length(power) > 1)
             || is.infinite((power)) || power < 0.0) {
    stop("Power provided is not a finite double  - distribution 
         cannot be constructed")
  }
}

