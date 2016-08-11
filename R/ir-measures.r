# Functions to compute IR measures for pairwise loss for
# a single group
# Notes:
# * Inputs are passed as a 2-element (y,f) list, to
#   facilitate the 'by' iteration
# * Return the respective metric, or a negative value if
#   it is undefined for the given group
# * For simplicity, we have no special handling for ties;
#   instead, we break ties randomly. This is slightly
#   inaccurate for individual groups, but should have
#   a small effect on the overall measure.


# Area under ROC curve = ratio of correctly ranking pairs


#' Compute Information Retrieval measures.
#' 
#' Functions to compute Information Retrieval measures for pairwise loss for a
#' single group. The function returns the respective metric, or a negative
#' value if it is undefined for the given group.
#' 
#' For simplicity, we have no special handling for ties; instead, we break ties
#' randomly. This is slightly inaccurate for individual groups, but should have
#' only a small effect on the overall measure.
#' 
#' \code{gbm_conc} computes the concordance index: Fraction of all pairs (i,j)
#' with i<j, x[i] != x[j], such that x[j] < x[i]
#' 
#' If \code{obs} is binary, then \code{gbm_roc_area(obs, pred) =
#' gbm.conc(obs[order(-pred)])}.
#' 
#' \code{gbm_conc} is more general as it allows non-binary targets, but is
#' significantly slower.
#' 
#' @aliases gbm_roc_area gbm_conc
#' @param obs Observed value
#' @param pred Predicted value
#' @return The requested performance measure.
#' @author Stefan Schroedl
#' @seealso \code{\link{gbm}}
#' @references C. Burges (2010). "From RankNet to LambdaRank to LambdaMART: An
#' Overview", Microsoft Research Technical Report MSR-TR-2010-82.
#' @keywords models
#' @examples
#' 
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--	or do  help(data=index)  for the standard data sets.
#' 
#' @export
gbm_roc_area <- function(obs, pred) {
   n1 <- sum(obs)
   n <- length(obs)
   if (n==n1) { return(1) }
   # Fraction of concordant pairs
   # = sum_{pos}(rank-1) / #pairs with different labels
   # #pairs = n1 * (n-n1)
   return ((mean(rank(pred)[obs > 0]) - (n1 + 1)/2)/(n - n1))
}

# Concordance Index:
# Fraction of all pairs (i,j) with i<j, x[i] != x[j], such that x[j] < x[i]
# Invariant: if obs is binary, then
#      gbm_roc_area(obs, pred) = gbm_conc(obs[order(-pred)])
# gbm_conc is more general as it allows non-binary targets,
# but is significantly slower
gbm_conc <- function(x) {
   lx <- length(x)
   return (sum(mapply(function(r) { sum(x[(r+1):lx]<x[r]) }, 1:(lx-1))))
}


