#' Adjust scales of prediction
#' 
#' S3 method (used internally within predict.GBMFit) that adjusts the scales of the predictor variable.
#' 
#' @usage adjust_pred_scale(gbm_predictions, distribution_obj)
#' 
#' @param gbm_predictions a matrix containing predictions generated from the C++ function \code{\link{gbm_pred}}.
#' 
#' @param distribution_obj a GBMDist object used in the fitting a gbm model and generating the gbm_predictions.
#' 
#' @return gbm_predictions scaled to the response of the distribution
#' 

adjust_pred_scale <- function(gbm_predictions, distribution_obj) {
  UseMethod("adjust_pred_scale", distribution_obj)
}

adjust_pred_scale.AdaBoostGBMDist <- function(gbm_predictions, distribution_obj) {
  gbm_predictions <- 1 / (1 + exp(-2*gbm_predictions))
  return(gbm_predictions)
}

adjust_pred_scale.BernoulliGBMDist <- function(gbm_predictions, distribution_obj) {
  gbm_predictions <- 1/(1+exp(-gbm_predictions))
  return(gbm_predictions)
}

adjust_pred_scale.CoxPHGBMDist <- function(gbm_predictions, distribution_obj) {
  return(gbm_predictions)
}

adjust_pred_scale.GammaGBMDist <- function(gbm_predictions, distribution_obj) {
  gbm_prediction <- exp(gbm_predictions)
  return(gbm_predictions)
}

adjust_pred_scale.GaussianGBMDist <- function(gbm_predictions, distribution_obj) {
  return(gbm_predictions)
}

adjust_pred_scale.HuberizedGBMDist <- function(gbm_predictions, distribution_obj) {
  return(gbm_predictions)
}

adjust_pred_scale.LaplaceGBMDist <- function(gbm_predictions, distribution_obj) {
  return(gbm_predictions)
}

adjust_pred_scale.PairwiseGBMDist <- function(gbm_predictions, distribution_obj) {
  gbm_predictions <- 1/(1+exp(-gbm_predictions))
  return(gbm_predictions)
}

adjust_pred_scale.PoissonGBMDist <- function(gbm_predictions, distribution_obj) {
  gbm_predictions <- exp(gbm_predictions)
  return(gbm_predictions)
}

adjust_pred_scale.QuantileGBMDist <- function(gbm_predictions, distribution_obj) {
  return(gbm_predictions)
}

adjust_pred_scale.TDistGBMDist <- function(gbm_predictions, distribution_obj) {
  return(gbm_predictions)
}

adjust_pred_scale.TweedieGBMDist <- function(gbm_predictions, distribution_obj) {
  gbm_prediction <- exp(gbm_predictions)
  return(gbm_predictions)
}