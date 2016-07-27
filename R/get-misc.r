#' Get Misc Parameters
#' 
#' S3 method to get the "misc" objects from GBM distributions
#' 
#' @usage get_misc(distribution_obj)
#' 
#' @param distribution_obj a GBMDist object
#' 
#' @return a list of additional parameters specific to the distribution
#' 

get_misc <- function(distribution_obj) {
  check_if_gbm_dist(distribution_obj)
  UseMethod("get_misc", distribution_obj)
}

get_misc.default <- function(distribution_obj) {
  stop("Distribution not recognised - can't get misc")
}

get_misc.AdaBoostGBMDist <- function(distribution_obj) {
  return(as.list(NA))
}

get_misc.BernoulliGBMDist <- function(distribution_obj) {
  return(as.list(NA))
}

get_misc.CoxPHGBMDist <- function(distribution_obj) {
  return(list(ties=distribution_obj$ties))
}

get_misc.GammaGBMDist <- function(distribution_obj) {
  return(as.list(NA))
}

get_misc.GaussianGBMDist <- function(distribution_obj) {
  return(as.list(NA))
}

get_misc.HuberizedGBMDist <- function(distribution_obj) {
  return(as.list(NA))
}

get_misc.LaplaceGBMDist <- function(distribution_obj) {
  return(as.list(NA))
}

get_misc.PairwiseGBMDist <- function(distribution_obj) {
  return(list("GroupsAndRanks"=c(distribution_obj$group, distribution_obj$max_rank)))
}

get_misc.PoissonGBMDist <- function(distribution_obj) {
  return(as.list(NA))
}

get_misc.QuantileGBMDist <- function(distribution_obj) {
  return(list(alpha=distribution_obj$alpha))
}

get_misc.TDistGBMDist <- function(distribution_obj) {
  return(list(df=distribution_obj$df))
}

get_misc.TweedieGBMDist <- function(distribution_obj) {
  return(list(power=distribution_obj$power))
}