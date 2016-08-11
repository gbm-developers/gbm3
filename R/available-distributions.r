#' Available Distributions
#' 
#' Returns a vector of the distributions implemented in the package.
#' 
#' @usage available_distributions()
#' 
#' @return a vector containing the names of the available distributions in the package.
#' 
#' @author James Hickey 
#' 
#' @export 

available_distributions <- function() {
  distributions <- c("AdaBoost", "Bernoulli", "CoxPH", "Gamma", "Gaussian",
                     "Huberized", "Laplace", "Pairwise", "Poisson", "Quantile",
                     "TDist", "Tweedie")
  return(distributions)
}