# Verify Weights Group Consistency
# Function that updates the weights in gbm_data 
# so they are consistent with the groupings.



weight_group_consistency <- function(gbm_data, distribution_obj) {
  check_if_gbm_data(gbm_data)
  check_if_gbm_dist(distribution_obj)
  
  if(distribution_obj$name == "Pairwise") {
    
    if(is.null(distribution_obj$group_name)) {
      stop("Pairwise distribution object has no group defined.")
    }
    
    # Check that weights are constant across groups
    if ((!missing(gbm_data$weights)) && (!is.null(gbm_data$weights)))
    {
      w.min <- tapply(gbm_data$weights, INDEX=distribution_obj$group_name, FUN=min)
      w.max <- tapply(gbm_data$weights, INDEX=distribution_obj$group_name, FUN=max)
      
      if (any(w.min != w.max))
      {
        stop("For distribution 'pairwise', all instances for the same group must have the same weight")
      }
      
      # Normalize across groups
      gbm_data$weights <- gbm_data$weights * length(w.min) / sum(w.min)
    }
  }
  return(gbm_data)
}
