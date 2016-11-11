# Validate data
# 
# Validates that the data is consistent with the distribution selected.
# 
# @usage validata_gbm_data(gbm_data_obj, distribution_obj)
# 
# @param gbm_data_obj A GBMData object created using the function "gbm_data".
# 
# @param distribution_obj A GBMDist object storing all the parameters which define the model
# distribution selected.
# 
# @author James Hickey
#
# @return a validated gbm_data_obj
# 

validate_gbm_data <- function(gbm_data_obj, distribution_obj) {
  check_if_gbm_data(gbm_data_obj)
  check_if_gbm_dist(distribution_obj)
  
  # Normalise Weights
  if(!(distribution_name(distribution_obj) %in% available_distributions())) {
    stop("Distribution not recognised - see available_distributions for
           supported distributions")
  }

  if(distribution_name(distribution_obj) != "Pairwise") {
    gbm_data_obj$weights <- gbm_data_obj$weights / mean(gbm_data_obj$weights)
  } 
  
  # Check offset
  gbm_data_obj$offset <- check_offset(gbm_data_obj$offset, gbm_data_obj$y, distribution_obj)
  
  # Check responses
  check_response_values(distribution_obj, gbm_data_obj$y)
  
  return(gbm_data_obj)
}
