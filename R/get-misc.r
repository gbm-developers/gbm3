#' Get Misc Parameters
#' 
#' S3 method to get the "misc" objects from GBM distributions
#' 
#' @usage get_misc(distribution_obj)
#' 
#' @param distribution_obj a GBMDist object
#' 
#' @return a list of parameters
#' 

get_misc <- function(distribution_obj) {
  check_if_gbm_dist(distribution_obj)
  UseMethod("get_misc", distribution_obj)
}

get_misc.AdaBoostGBMDist <- function(distribution_obj) {
  return(as.list(NA))
}

get_misc.BernoulliGBMDist <- function(distribution_obj) {
  return(as.list(NA))
}

get_misc.CoxPHGBMDist <- function(distribution_obj) {
  return(as.list(c(ties=distribution_obj$ties)))
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
  return(as.list(Na))
}

get_misc.PairwiseGBMDist <- function(distribution_obj) {
  return(list("GroupsAndRanks"=c(distribution_obj$group, distribution_obj$max_rank)))
}

get_misc.PoissonGBMDist <- function(distribution_obj) {
  return(as.list(NA))
}

get_misc.QuantileGBMDist <- function(distribution_obj) {
  return(c(alpha=distribution_obj$alpha))
}

get_misc.TDistGBMDist <- function(distribution_obj) {
  return(as.list(distribution_obj$df))
}

get_misc.TweedieGBMDist <- function(distribution_obj) {
  return(as.list(c(power=distribution_obj$power)))
}