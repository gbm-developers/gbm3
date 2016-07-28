#' Create distributions
#' 
#' Creation methods that take an empty distribution object
#' and builds on top of the skeleton.
#' 
#' @usage create_dist(empty_obj, ...)
#' 
#' @param empty_obj  An empty distribution object of the correct class.
#' 
#' @param ...  extra parameters used to define the distribution object,
#'  see \code{\link{gbm_dist}}.
#'  
#' @return an appropriated gbm distribution object
#' 

create_dist <- function(empty_obj, ...) {
  check_if_gbm_dist(empty_obj)
  UseMethod("create_dist", empty_obj)
}

create_dist.default <- function(empty_obj, ...) {
  stop("couldn't create distribution object as creation method for selected distribution not found.")
}

create_dist.AdaBoostGBMDist <- function(empty_obj, ...) {
  check_dist_params(empty_obj, ...)
  return(empty_obj)
}

create_dist.BernoulliGBMDist <- function(empty_obj, ...) {
  check_dist_params(empty_obj, ...)
  return(empty_obj)
}

create_dist.CoxPHGBMDist <- function(empty_obj, strata=NA, sorted=NA, ties="efron"
                                     , prior_node_coeff_var=1000, ...) {
  check_dist_params(empty_obj, strata, sorted, ties, prior_node_coeff_var, ...)
  if(!(ties %in% c("breslow", "efron"))) {
    warning("Ties method not recognised - defaulting to efron")
    ties <- "efron"
  }
  empty_obj$ties <- ties
  empty_obj$original_strata_id <- convert_strata(strata)
  empty_obj$sorted <- sorted
  empty_obj$prior_node_coeff_var <- prior_node_coeff_var
  empty_obj$reorder <- TRUE
  return(empty_obj)
}

create_dist.GammaGBMDist <- function(empty_obj, ...) {
  check_dist_params(empty_obj, ...)
  return(empty_obj)
}

create_dist.GaussianGBMDist <- function(empty_obj, ...) {
  check_dist_params(empty_obj, ...)
  return(empty_obj)
}

create_dist.HuberizedGBMDist <- function(empty_obj, ...) {
  check_dist_params(empty_obj, ...)
  return(empty_obj)
}

create_dist.LaplaceGBMDist <- function(empty_obj, ...) {
  check_dist_params(empty_obj, ...)
  return(empty_obj)
}

create_dist.PairwiseGBMDist <- function(empty_obj, group="query", metric="ndcg", max_rank=0,
                                        ...) {
  check_dist_params(empty_obj, group, metric, max_rank, ...)
  empty_obj$metric <- metric
  empty_obj$group <- group
  empty_obj$max_rank <- max_rank
  empty_obj$reorder <- TRUE
  return(empty_obj)
}

create_dist.PoissonGBMDist <- function(empty_obj, ...) {
  check_dist_params(empty_obj, ...)
  return(empty_obj)
}

create_dist.QuantileGBMDist <- function(empty_obj, alpha=0.25, ...) {
  check_dist_params(empty_obj, alpha, ...)
  empty_obj$alpha <- alpha 
  return(empty_obj)
}

create_dist.TDistGBMDist <- function(empty_obj, df=4, ...) {
  if(is.null(df)) df <- 4
  check_dist_params(empty_obj, df, ...)
  empty_obj$df <- df
  return(empty_obj)
}

create_dist.TweedieGBMDist <- function(empty_obj, power=1.5, ...) {
  if(is.null(power)) power <- 1.5
  check_dist_params(empty_obj, power, ...)
  empty_obj$power <- power
  return(empty_obj)
}
