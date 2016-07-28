#' Set Up Groups
#' 
#' Internal S3 method to set up groups - currently this only
#' applies to the pairwise distribution. 
#' 
#' @usage determine_groups(original_data, response, distribution_obj)
#' 
#' @param original_data the original data used on the initial call to \code{\link{gbm2}}.
#' 
#' @param response the response data provided on the initial call to \code{\link{gbm2}}.
#' 
#' @param distribution_obj a GBMDist object created using gbm_dist.
#' 
#' @return an updated distribution_object with the group, group_order and group_index updated
#' 
#' @export determine_groups

determine_groups <- function(original_data, response, distribution_obj) {
  check_if_gbm_dist(distribution_obj)
  UseMethod("determine_groups", distribution_obj) 
}

determine_groups.default <- function(original_data, response, distribution_obj) {
  return(distribution_obj)
}

determine_groups.PairwiseGBMDist <- function(original_data, response, distribution_obj) {
  if (is.null(distribution_obj$group)) {
    stop("For pairwise regression, the distribution parameter must be a list with a parameter 'group' for 
         the list of the column names indicating groups, for example group=c(\"date\",\"session\",\"category\",\"keywords\").")
  }
  
  # Check if grouping is specified in data
  i <- match(distribution_obj$group, colnames(original_data))
  if (any(is.na(i))) {
    stop("Group column does not occur in data: ", distribution_obj$group[is.na(i)])
  }
  
  # Construct group index
  distribution_obj$group_index <- factor(do.call(paste, c(original_data[,distribution_obj$group, drop=FALSE], sep=":")))

  # Shuffle groups, to remove bias when splitting into train/test set and/or CV folds
  perm_levels  <- levels(distribution_obj$group_index)[sample(1:nlevels(distribution_obj$group_index))]
  distribution_obj$group        <- factor(distribution_obj$group_index, levels=perm_levels)
  
  # The C++ function expects instances to be sorted by group and descending by target
  distribution_obj$group_order    <- order(distribution_obj$group, -response)
  distribution_obj$group  <- distribution_obj$group[distribution_obj$group_order]

  return(distribution_obj)
}