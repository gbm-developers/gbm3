# Function that updates the number of training rows and observations
# based off of groups
# @author James Hickey

update_num_train_groups <- function(train_params, dist_obj) {
  check_if_gbm_train_params(train_params)
  check_if_gbm_dist(dist_obj)
  UseMethod("update_num_train_groups", dist_obj)
}

update_num_train_groups.default <- function(train_params, dist_obj) {
  return(train_params)
}

update_num_train_groups.PairwiseGBMDist <- function(train_params, dist_obj) {
  # Check if more than one row per observation
  if(any(duplicated(train_params$id))) {
    stop("Pairwise does not currently support multiple rows per observation.")
  }
  
  # Calculate number of training groups
  num_train_groups <- max(1, round(train_params$train_fraction * nlevels(dist_obj$group)))
  
  # Get number of training rows - include those up to max group
  train_params$num_train_rows <- max(which(dist_obj$group == levels(dist_obj$group)[num_train_groups]))
  train_params$num_train <- train_params$num_train_rows
  
  return(train_params)
}