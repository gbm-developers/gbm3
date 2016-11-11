# Verify Weights Group Consistency
# Function that updates the weights in gbm_data 
# so they are consistent with the groupings.

weight_group_consistency <- function(distribution_obj, weights) {
    UseMethod("weight_group_consistency")
}

weight_group_consistency.default <- function(distribution_obj, weights) {
    stop("The impossible has happened")
}

weight_group_consistency.GBMDist <- function(distribution_obj, weights) {
    weights
}

weight_group_consistency.PairwiseGBMDist <-
    function(distribution_obj, weights) {
        if(is.null(distribution_obj$group_index)) {
            stop("Pairwise distribution object has no group defined.")
        }
        
        ## Check that weights are constant across groups
        if ((!missing(weights)) && (!is.null(weights))) {
            w_min <- tapply(weights,
                            INDEX=distribution_obj$group_index, FUN=min)
            w_max <- tapply(weights,
                            INDEX=distribution_obj$group_index, FUN=max)
      
            if (any(w_min != w_max)) {
                stop("For distribution 'pairwise', all instances for the same group must have the same weight")
            }
      
      # Normalize across groups
            weights <- weights * length(w_min) / sum(w_min)
        }

        weights
    }
