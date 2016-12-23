
#' Check, validate, and order gbm data
#' @inheritParams gbmt_fit
#' @return \code{GBMPrep} object
#' @export
gbmt_data <- function(x, y, distribution=gbm_dist("Gaussian"),
                       weights=rep(1, nrow(x)), offset=rep(0, nrow(x)),
                       train_params=training_params(num_trees=100,
                                                    interaction_depth=3,
                                                    min_num_obs_in_node=10,
                                                    shrinkage=0.001, bag_fraction=0.5,
                                                    id=seq_len(nrow(x)), num_train=round(0.5 * nrow(x)),
                                                    num_features=ncol(x)), response_name="y",
                       var_monotone=NULL,
                       var_names=NULL, keep_gbm_data=FALSE, cv_folds=1,
                       cv_class_stratify=FALSE, fold_id=NULL,
                       par_details=getOption('gbm.parallel'),
                       is_verbose=FALSE){
  
  # Check num_features makes sense
  if(train_params$num_features > ncol(x)) {
    warning("Number of features exceeds the number of predictor variables - 
            setting number of features to number of predictors")
    train_params$num_features <- ncol(x)
  }
  
  # Create gbm_data_obj
  gbm_data_obj <- gbm_data(x, y, weights, offset)
  
  # Set-up variable containers
  variables <- var_container(gbm_data_obj, var_monotone, var_names)
  
  # Create strata
  distribution <- create_strata(gbm_data_obj, train_params, distribution)
  
  # Process data obj and validate
  gbm_data_obj <- convert_factors(gbm_data_obj)
  gbm_data_obj <- validate_gbm_data(gbm_data_obj, distribution)
  
  # Order the data
  gbm_data_obj <- order_data(gbm_data_obj, distribution, train_params)
  
  # Order ids according to group_order if necessary
  if(!is.null(distribution$group_order)) {
    train_params$id <- train_params$id[distribution$group_order]
  } 
  train_params$id <- train_params$id[train_params$id_order]
  
  # Get CV groups
  cv_groups <- create_cv_groups(gbm_data_obj, distribution,
                                train_params, cv_folds,
                                cv_class_stratify, fold_id)
  
  gbmt_prep <- list(gbm_data_obj=gbm_data_obj,
              distribution=distribution,
              train_params=train_params,
              variables=variables,
              cv_folds=cv_folds,
              cv_groups=cv_groups,
              par_details=par_details,
              is_verbose=is_verbose,
              keep_gbm_data=keep_gbm_data,
              fold_id=fold_id,
              response_name=response_name)
  class(gbmt_prep) <- "GBMPrep"
  return(gbmt_prep)
}

#' Fit a gbm model to a gbm prep obj
#' @param object GBMPrep object from \code{gbmt_data}
#' @export
gbmt_fit_ <- function(object){
  # Create fitted object
  gbm_fit_obj <- gbm_cross_val(object$gbm_data_obj, object$distribution,
                               object$train_params, object$variables,
                               object$cv_folds, object$cv_groups,
                               object$par_details, object$is_verbose)
  
  # Keep original data
  if(object$keep_gbm_data) {
    gbm_fit_obj$gbm_data_obj <- object$gbm_data_obj
  }
  # Reorder if necessary
  gbm_fit_obj <- reorder_fit(gbm_fit_obj, object$distribution)
  gbm_fit_obj$fold_ids <- object$fold_id
  gbm_fit_obj$par_details <- object$par_details
  gbm_fit_obj$response_name <- object$response_name
  
  return(gbm_fit_obj)
} 
