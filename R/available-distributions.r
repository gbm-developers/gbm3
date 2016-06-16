#' Available Distributions
#' 
#' Returns a vector of the distributions implemented in the package.
#' 
#' @usage available_distributions()
#' 
#' @return a vector o
#' 
#' @export available_distributions

available_distributions <- function() {
  distributions <- c("AdaBoost", "Bernoulli", "CoxPH", "Gamma", "Gaussian",
                     "Huberized", "Laplace", "Pairwise", "Poisson", "Quantile",
                     "TDist", "Tweedie")
  return(distributions)
}