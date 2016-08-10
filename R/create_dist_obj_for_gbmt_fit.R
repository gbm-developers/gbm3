# Extract parameters from gbm to create \code{GBMDist}
#
# Internal function which creates the appropriate \code{GBMDist} from 
# the arguments of a call to \code{gbm} or \code{gbm.fit}.
#
# @usage create_dist_obj_for_gbmt_fit(distribution, tied.times.method="efron", strata=NA, prior.node.coeff.var=1000)
#
# @param distribution the list of parameters defining the distribution for gbmt_fit - see \code{gbm}
#
# @param tied.times.method This is an optional string used with
# \code{CoxPH} distribution specifying what method to employ when dealing with tied times. 
# Currently only "efron" and "breslow" are available; the default value is "efron". 
# Setting the string to any other value reverts the method to the original CoxPH model
# implementation where ties are not explicitly dealt with.
#
# @param strata Optional vector of integers (or factors) only used with the \code{CoxPH} distributions.
# Each integer in this vector represents the stratum the corresponding row in the data belongs to,
# e. g. if the 10th element is 3 then the 10th data row belongs to the 3rd strata.
#
# @param prior.node.coeff.var Optional double only used with the \code{CoxPH} distribution.
# It is a prior on the coefficient of variation associated with the hazard 
# rate assigned to each terminal node when fitting a tree.Increasing its value emphasises 
# the importance of the training data in the node when assigning a prediction to said node.
#
# @return \code{GBMDist} object
# 
# @author James Hickey

create_dist_obj_for_gbmt_fit <- function(distribution, tied.times.method="efron", strata=NA, prior.node.coeff.var=1000) {
  # Check inputs
  if(is.null(distribution$name)) stop("distribution parameter name not defined.")
  if(!(tolower(distribution$name) %in% tolower(available_distributions()))) {
    stop("The distribution ", distribution$name, " is not available in the gbm package.")
  }
  
  # Create distribution
  if(tolower(distribution$name) == "coxph") {
    dist <- create_coxph_for_gbmt_fit(distribution, tied.times.method, strata, prior.node.coeff.var)
  } else if(tolower(distribution$name) == "pairwise") {
    dist <- create_pairwise_for_gbmt_fit(distribution)
  } else if(tolower(distribution$name) == "quantile") {
    dist <- create_quantile_for_gbmt_fit(distribution)
  } else if(tolower(distribution$name) == "tdist") {
    dist <- create_tdist_for_gbmt_fit(distribution)
  } else if(tolower(distribution$name) == "tweedie") {
    dist <- create_tweedie_for_gbmt_fit(distribution)
  } else {
    # Find distribution in list of available distributions
    match_dist <- match(tolower(distribution$name), tolower(available_distributions()))
    dist <- create_dist(empty_distribution(available_distributions()[match_dist]))
  }
  
  return(dist)
}

create_coxph_for_gbmt_fit <- function(distribution, tied.times.method, strata, prior.node.coeff.var) {
  return(create_dist(empty_distribution("CoxPH"), strata, sorted=NA, 
                     ties=tied.times.method, prior_node_coeff_var=prior.node.coeff.var))
}

create_pairwise_for_gbmt_fit <- function(distribution) {
  if(is.null(distribution$group)) distribution$group <- "query"
  if(is.null(distribution$metric)) distribution$metric <- "ndcg"
  if(is.null(distribution$max.rank)) distribution$max.rank <- 0
  return(create_dist(empty_distribution("Pairwise"), 
                     distribution$group, distribution$metric, distribution$max.rank, distribution$group_index))
}

create_quantile_for_gbmt_fit <- function(distribution) {
  return(create_dist(empty_distribution("Quantile"), distribution$alpha))
}

create_tdist_for_gbmt_fit <- function(distribution) {
  return(create_dist(empty_distribution("TDist"), distribution$df))
}

create_tweedie_for_gbmt_fit <- function(distribution) {
  return(create_dist(empty_distribution("Tweedie"), distribution$power))
}