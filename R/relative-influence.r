#' Methods for estimating relative influence
#' 
#' Helper functions for computing the relative influence of each variable in
#' the \code{GBMFit} object
#' 
#' These functions offer the different
#' methods for computing the relative influence in \code{\link{summary.GBMFit}}.
#'
#'@usage relative_influence(gbm_fit_obj, rescale, sort_it)
#'
#' @param gbm_fit_obj a \code{GBMFit} object created from an initial call to
#' \code{\link{gbm2}}.
#' 
#' @param num_trees the number of trees to use for computations. If not provided,
#' the function will guess: if a test set was used in fitting, the number of
#' trees resulting in lowest test set error will be used; otherwise, if
#' cross-validation was performed, the number of trees resulting in lowest
#' cross-validation error will be used; otherwise, all trees will be used.
#' 
#' @param rescale  whether or not the result should be scaled. Defaults to
#' \code{FALSE}.
#' 
#' @param sort_it  whether or not the results should be (reverse) sorted.
#' Defaults to \code{FALSE}.
#' 
#' @return By default, returns an unprocessed vector of estimated relative
#' influences. If the \code{rescale} and \code{sort} arguments are used,
#' returns a processed version of the same.
#' 
#' @author Greg Ridgeway \email{gregridgeway@@gmail.com}
#' @details \code{\link{relative.influence}} is the same as that
#' described in Friedman (2001).
#' \code{\link{permutation.test.gbm}} randomly permutes each
#' predictor variable at a time and computes the associated reduction in
#' predictive performance. This is similar to the variable importance measures
#' Breiman uses for random forests, but \code{gbm2} currently computes using the
#' entire training dataset (not the out-of-bag observations).

#' @seealso \code{\link{summary.gbm}}
#' @references J.H. Friedman (2001). "Greedy Function Approximation: A Gradient
#' Boosting Machine," Annals of Statistics 29(5):1189-1232.
#' 
#' L. Breiman (2001).
#' \href{http://oz.berkeley.edu/users/breiman/randomforest2001.pdf}{Random
#' Forests}.
#' @keywords hplot
#' @export relative_influence
#'

relative_influence <- function(gbm_fit_obj, num_trees, rescale = FALSE, sort_it = FALSE){
  # Initial checks
  check_if_gbm_fit(gbm_fit_obj)
  if(!is.logical(rescale) || (length(rescale) > 1))
    stop("rescale argument must be a logical")
  if(!is.logical(sort_it) || (length(sort_it) > 1))
    stop("sort_it must be a logical")
  
  # Fill in missing values
  if( missing( num_trees ) ){
    if ( gbm_fit_obj$params$train_fraction < 1 ){
      num_trees <- gbm_perf(gbm_fit_obj, method="test", plot_it=FALSE )
    }
    else if ( !is.null( gbm_fit_obj$cv_error ) ){
      num_trees <- gbm_perf( gbm_fit_obj, method="cv", plot_it = FALSE )
    }
    else{
      num_trees <- gbm_fit_obj$params$num_trees
    }
    message(paste( "num_trees not given. Using", num_trees, "trees.\n" ))
    
  }

  # Create relative influence for every variable
  rel_inf_verbose <- unlist(lapply(gbm_fit_obj$trees[seq_len(num_trees)], get_rel_inf_of_vars))
  
  # Sum across trees and remove unused variables (names are "-1")
  rel_inf_compact <- unlist(lapply(split(rel_inf_verbose, names(rel_inf_verbose)), sum))
  rel_inf_compact <- rel_inf_compact[names(rel_inf_compact) != "-1"]
  
  # rel_inf_compact excludes those variable that never entered the model
  # insert 0's for the excluded variables
  rel_inf <- rep(0, length(gbm_fit_obj$variables$var_names))
  i <- as.numeric(names(rel_inf_compact))+1
  rel_inf[i] <- rel_inf_compact
  names(rel_inf) <- gbm_fit_obj$variables$var_names
  
  # Rescale and sort
  if (rescale) rel_inf <- rel_inf / max(rel_inf)
  if (sort_it) rel_inf <- rev(sort(rel_inf))
  
  return(rel_inf)
}

#### Helper function ####
get_rel_inf_of_vars <- function(obj) {
  lapply(split(obj[[6]], obj[[1]]), sum) # 6 - Improvement, 1 - var name
}