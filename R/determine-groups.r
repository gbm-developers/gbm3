#' Set Up Groups
#' 
#' Function to set up groups - currently this only
#' applies to the pairwise distribution
#' 
#' @usage  determine_groups(col_names, gbm_data, distribution_obj)
#' 
#' @param col_names vector of variable names used in fittig - required for Pairwise
#' @param gbm_data a GBMData object created using 
#' 
#' @param distribution_obj a GBMDist object created using gbm_dist
#' 
#' @return an updated distribution_object 
#' 

determine_groups <- function(col_names, gbm_data, distribution_obj) {
  # Check names
  if(!is.atomic(col_names) || any(col_names != as.character(col_names))
     || is.null(col_names)) {
    stop("Names of data must be a vector of strings.")
  }
  
  check_if_gbm_data(gbm_data)
  check_if_gbm_dist(distribution_obj)
  if(distribution_obj$reorder) {
   UseMethod("determine_groups", distribution_obj) 
  }
  
  return(distribution_obj)
}

determine_groups.CoxPHGBMDist <- function(col_names, gbm_data, distribution_obj) {
  return(distribution_obj)
}

determine_groups.PairwiseGBMDist <- function(col_names, gbm_data, distribution_obj) {
  if (is.null(distribution_obj$group))
  {
    stop("For pairwise regression, the distribution parameter must be a list with a parameter 'group' for 
         the list of the column names indicating groups, for example group=c(\"date\",\"session\",\"category\",\"keywords\").")
  }
  
  # Check if grouping is specified in data
  i <- match(distribution_obj$group, col_names)
  if (any(is.na(i)))
  {
    stop("Group column does not occur in data: ", distribution_obj$group[is.na(i)])
  }
  
  # Construct group index
  distribution_obj$group <- factor(do.call(paste, c(gbm_data[,distribution_obj$group, drop=FALSE], sep=":")))

  # Shuffle groups, to remove bias when splitting into train/test set and/or CV folds
  perm.levels  <- levels(distribution_obj$group)[sample(1:nlevels(distribution_obj$group))]
  distribution_obj$group        <- factor(distribution_obj$group, levels=perm.levels)
  
  # The C++ function expects instances to be sorted by group and descending by target
  distribution_obj$group_order    <- order(distribution_obj$group, -gbm_data$y)
  distribution_obj$group  <- distribution_obj$group[distribution_obj$group_order]

  return(distribution_obj)
}