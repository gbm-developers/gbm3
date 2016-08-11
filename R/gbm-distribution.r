#' GBM Distribution
#' 
#' Generates distribution object for gbmt.
#' 
#' @usage gbm_dist(name="Gaussian", ...)
#' 
#' @param name The name (a string) of the distribution to be initialized and used
#' in fitting a gradient boosted model via gbmt.  The current distributions
#' available can be viewed using the available-distributions function.  If no 
#' distribution is specified this function constructs a Gaussian distribution by
#' default.
#' 
#' @param ...  Extra named parameters required for initializing certain distributions.  
#' If t-distribution is selected, an additional parameter (\code{df}) specifying the number of 
#' degrees of freedom can be given.  The default degrees of freedom is set to four.  
#' 
#' If quantile is seleceted then the quantile to estimate may be specified using the
#' named parameter \code{alpha}.  The default quantile to estimate is 0.25.
#' 
#' If the tweedie distribution is selected the power-law specifying the distribution
#' may be set via the named parameter: \code{power}.  This parameter defaults to unity.
#' 
#' If a Cox Partial Hazards model is selected a number of additional parameters are
#' required, these are:
#' 
#' \describe{
#'  \item{\code{strata}}{A vector of integers (or factors) specifying which strata each data-row belongs to, 
#'   if none is specified it is assumed all training data is in the same stratum.}
#'   \item{\code{ties}}{String specifying the method to be used when dealing with tied
#'   event times.  Currently only "breslow" and "efron" are available, with the latter
#'   being the default.}
#'   \item{\code{prior_node_coeff_var}}{It is a prior on the coefficient of variation associated with the hazard 
#' rate assigned to each terminal node when fitting a tree.Increasing its value emphasises 
#' the importance of the training data in the node when assigning a prediction to said node. This defaults to 1000.}
#' }
#' 
#' Finally, if the pairwise distribution is selected a number of parameters also need to be
#' specified. These parameters are \code{group}, \code{metric} and \code{max_rank}.
#' The first is a character vector  with the column names of data that jointly indicate the group 
#' an instance bleongs to (typically a query in Information Retrieval).  For
#' training, only pairs of instances from the same group and with different target
#' labels may be considered. \code{metric} is
#' the IR measure to use, one of
#' \describe{
#'   \item{list("conc")}{Fraction of concordant pairs; for binary labels,
#'    this is equivalent to the Area under the ROC Curve}
#'   \item{:}{Fraction of concordant pairs; for binary labels, this
#' is equivalent to the Area under the ROC Curve}
#'   \item{list("mrr")}{Mean reciprocal rank of the highest-ranked positive instance}
#'   \item{:}{Mean reciprocal rank of the highest-ranked positive instance}
#'   \item{list("map")}{Mean average precision, a generalization of \code{mrr}
#'   to multiple positive instances}
#'   \item{:}{Mean average precision, a generalization of \code{mrr} to multiple
#'    positive instances}
#'   \item{list("ndcg:")}{Normalized discounted cumulative gain.
#' The score is the weighted sum (DCG) of the user-supplied target values, weighted by
#' log(rank+1), and normalized to the maximum achievable value. This is the
#' default if the user did not specify a metric.}
#' }
#' 
#' \code{ndcg} and \code{conc} allow arbitrary target values, while binary
#' targets {0,1} are expected for \code{map} and \code{mrr}. For \code{ndcg}
#' and \code{mrr}, a cut-off can be chosen using a positive integer parameter
#' \code{max_rank}. If left unspecified, all ranks are taken into account.
#' 
#' @author James Hickey
#' 
#' @return returns a \code{GBMDist} object.
#' 
#' @export
#' 

gbm_dist <- function(name="Gaussian", ...) {
  # Check if name is available
  if(!(name %in% available_distributions())) {
    stop("The distribution ", name, " is not available in the gbm package.")
  }
  
  # Get a skeleton distribution object
  bare_distribution <- empty_distribution(name)
  
  # Create object and return
  return(create_dist(bare_distribution, ...))
}

