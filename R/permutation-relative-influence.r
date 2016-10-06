#' Relative influence via permutation
#'
#' This function offers a method for computing the relative influence in
#' \code{\link{summary.GBMFit}}, and is not intended to be called directly.
#'
#' Calculates the relative influence of predictors via random
#' permutation of each predictor one at a time and calculating the
#' associated reduction in predictive performance.  This experimental
#' measure is similar to the variable importance measures Breiman uses
#' for random forests, but \code{\link{gbmt}} currently computes using
#' the entire training dataset (not the out-of-bag observations).
#' 
#' @param gbm_fit_obj a \code{GBMFit} object from an initial call to
#' \code{\link{gbmt}}.
#' 
#' @param num_trees the number of trees to use for computations. If
#' not provided, the function will guess: if a test set was used in
#' fitting, the number of trees resulting in lowest test set error
#' will be used; otherwise, if cross-validation was performed, the
#' number of trees resulting in lowest cross-validation error will be
#' used; otherwise, all trees will be used.
#' 
#' @param rescale whether or not the result should be scaled. Defaults
#' to \code{FALSE}.
#' 
#' @param sort_it whether or not the results should be (reverse)
#' sorted.  Defaults to \code{FALSE}.
#' 
#' @return By default, returns an unprocessed vector of estimated
#' relative influences. If the \code{rescale} and \code{sort}
#' arguments are used, returns a processed version of the same.
#' 
#' @author Greg Ridgeway \email{gregridgeway@@gmail.com}
#'
#' @seealso \code{\link{summary.GBMFit}}
#' \href{https://www.stat.berkeley.edu/~breiman/randomforest2001.pdf}{Random
#' Forests}.
#' @keywords hplot
#' @export
permutation_relative_influence <- function(gbm_fit_obj,
                                           num_trees,
                                           rescale=FALSE,
                                           sort_it=FALSE) {
  # Checks initial inputs
  check_if_gbm_fit(gbm_fit_obj)
  check_if_natural_number(num_trees)
  if(num_trees > length(gbm_fit_obj$trees))
    stop("num_trees exceeds maximum")
  if(is.na(rescale) || !is.logical(rescale) || (length(rescale) > 1))
    stop("rescale argument must be a logical")
  if(is.na(sort_it) || !is.logical(sort_it) || (length(sort_it) > 1))
    stop("sort_it must be a logical")
  if(is.null(gbm_fit_obj$gbm_data_obj))
    stop("Model was fit with keep_data=FALSE: permutation_relative_influence has not been implemented for that case.")
  
  # Get variables used in the model
  variables_indices <- sort(unique(unlist(lapply(gbm_fit_obj$trees[seq_len(num_trees)],
                                      function(x){unique(x[[1]])}))))
  
  # Remove unused variables (those = "-1") and adjust indices from C convention 
  variables_indices <- variables_indices[variables_indices != -1] + 1
  
  # Set up rel_inf results
  rel_inf <- rep(0, length(gbm_fit_obj$variables$var_names))
  
  # Extract data for loss calculation
  y            <- gbm_fit_obj$gbm_data_obj$y
  offset       <- gbm_fit_obj$gbm_data_obj$offset
  weights      <- gbm_fit_obj$gbm_data_obj$weights
  x            <- matrix(gbm_fit_obj$gbm_data_obj$x, ncol=length(gbm_fit_obj$variables$var_names))
  gbm_fit_obj$Terms <- NULL # this makes predict.GBMFit take x as it is
    
  # Shuffle rows - indices
  shuffled_rows <- sample(1:nrow(x))
  for(i in seq_len(length(variables_indices))) {
    # Shuffle  predictor variables
    x[ ,variables_indices[i]]  <- x[shuffled_rows, variables_indices[i]]
    new_preds <- predict(gbm_fit_obj,
                         newdata=as.data.frame(x),
                         n.trees=num_trees)
    
    # Calculate loss with shuffled data
    rel_inf[variables_indices[i]] <- loss(y, new_preds, weights, offset,
                                          gbm_fit_obj$distribution,
                                          baseline=rep(gbm_fit_obj$train.error[num_trees], length(weights)))
    
    # Unshuffle variable - only permute variables one at a time
    x[shuffled_rows, variables_indices[i]] <- x[ ,variables_indices[i]]
  }
  
  # Scale and sort relative influence
  if (rescale) rel_inf <- rel_inf / max(rel_inf)
  if (sort_it) rel_inf <- rev(sort(rel_inf))
  
  return(rel_inf)
}
