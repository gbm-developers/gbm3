#' GBM Cross Validation
#' 
#' Create fitted object - Cross validation is performed if \code{cv_folds > 1}, otherwise only the full model
#' is fitted. This function is called by \code{\link{gbm2}} to perform the model fitting.
#' 
#' @usage gbm_cross_validation(gbm_data_obj, gbm_dist_obj, train_params, 
#' var_container, cv_folds, cv_groups, is_verbose)
#' 
#' @param gbm_data_obj a GBMData object containing all of the data used to fit a gbm model. 
#' 
#' @param gbm_dist_obj a GBMDist object specifying the distribution and any additional parameters needed.
#' 
#' @param train_params a GBMTrainParams object containing generic parameters defining how the model should be
#' trained.
#' 
#' @param var_container a GBMVarCont object which defines the properties of the predictor variables in the data.
#' 
#' @param cv_folds a positive integer specifying the number of folds to be used in cross-validation of the gbm
#' fit. If cv_folds > 1 then cross-validation is performed; the default of cv_folds is 1.
#' 
#' @param cv_groups vector of integers specifying which row of data belongs to which cv_folds.
#' 
#' @param par_details Details of the parallelization to use in the
#'     core algorithm.
#'     
#' @param is_verbose if TRUE, will print out progress and performance of the fitting.
#' 
#' @return a \code{GBMFit} object which contains appropriate CV info if requested.
#' 

gbm_cross_val <- function(gbm_data_obj, gbm_dist_obj, train_params, var_container, 
                          cv_folds, cv_groups, par_details, is_verbose) {
  # Create output list
  gbm_results <- list() 
  
  # Full model fit
  if(is_verbose) message("Fitting Final Model \n")
  gbm_results[[length(gbm_results)+1]] <- gbm_call(gbm_data_obj, gbm_dist_obj, train_params,
                                                  var_container, par_details, is_verbose)

  # Check if only need to fit full model
  if(cv_folds == 1) {
    gbm_results[[1]]$cv_folds <- 1
    return(gbm_results[[1]])
  }

  # Loop over folds
  for(fold_num in seq_len(cv_folds)) {
    if(is_verbose) message("CV:", fold_num, "\n")
    
    # Extract observations in cv fold
    gbm_object_list <- extract_obs_in_fold(gbm_data_obj, gbm_dist_obj, train_params, cv_groups, fold_num)
    
    # Fit to fold
    gbm_results[[length(gbm_results)+1]] <- gbm_call(gbm_object_list$data, gbm_object_list$dist, 
                                                    gbm_object_list$params, var_container, par_details, is_verbose)
  }
  
  # If have multiple folds then total results object is a different class
  class(gbm_results) <- "GBMCVFit"  
  
  # Calculate errors
  cv_errors <- gbm_cv_errors(gbm_results, cv_folds, cv_groups)
  
  # Best number of trees
  best_iter_cv <- which.min(cv_errors)
  
  # Calculate predictions
  predictions <- predict(gbm_results, gbm_data_obj, cv_folds, cv_groups, best_iter_cv)
  
  # Extract relevant parts - all data model
  gbm_results <- gbm_results[[1]]
  gbm_results$cv_folds <- cv_folds
  gbm_results$cv_error <- cv_errors
  gbm_results$cv_fitted <- predictions
  
  return(gbm_results)
}