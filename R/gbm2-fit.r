#' GBM2 fit
#' 
#' Fits a generalized boosting model.  This is for "power" users who have a large number of 
#' variables who wish to avoid calling \code{model.frame} which can be slow in this instance.
#' 
#' @usage gbm2.fit(x, y, distribution=gbm_dist("Gaussian", ...), weights=rep(1, nrow(x)), offset=rep(0, nrow(x))
#' train_params=training_params(num_trees=100, interaction_depth=1, min_num_obs_in_node=10, 
#' shrinkage=0.001, bag_fraction=0.5, id=seq_len(nrow(x)), num_train=round(0.5 * nrow(x)), num_features=ncol(x)), 
#' var_monotone=NULL, var_names=NULL, keep_gbm_data=FALSE, cv_folds=1, cv_class_stratify=FALSE, fold_id=NULL, is_verbose=FALSE)
#' 
#' @param x a data frame or data matrix containing the predictor variables. 
#' 
#' @param y is a matrix of outcomes. Excluding CoxPH this matrix of outcomes collapses to a vector; in the case
#' of CoxPH it is a survival object where the event times fill the first one (or two columns) and the statys fills the final
#' column.  The length of the 1st dimension of y must match the number of rows in x.
#' 
#' @param distribution a GBMDist object specifying the distribution and any additional parameters needed.
#' 
#' @param weights optional vector of weights used in the fitting process.  These weights must be positive but 
#' need not be normalized. By default they are set to 1 for each data row.
#' 
#' @param offset  optional vector specifying the model offset; must be positive.  This defaults to a vector of 0's, the length
#' of which is equal to the rows of x.
#' 
#' @param train_params  a GBMTrainParams object which specifies the parameters used in growing decision trees.
#' 
#' @param var_monotone optional vector, the same length as the number of predictors, indicating the relationship
#' each variable has with the outcome.  It have a monotone increasing (+1) or decreasing (-1) or an arbitrary relationship.
#' 
#' @param var_names a vector of strings of containing the names of the predictor variables.
#' 
#' @param keep_gbm_data a bool specifying whether or not the gbm_data object created in this method should be
#' stored in the results.
#' 
#' @param cv_folds a positive integer specifying the number of folds to be used in cross-validation of the gbm
#' fit. If cv_folds > 1 then cross-validation is performed; the default of cv_folds is 1.
#' 
#' @param cv_class_stratify a bool specifying whether or not to stratify via response outcome. Currently only 
#' applies to "Bernoulli" distribution and defaults to false.
#' 
#' @param fold_id An optional vector of values identifying what fold
#' each observation is in. If supplied, cv_folds can be missing. Note:
#' Multiple rows of the same observation must have the same fold_id.
#' 
#' @param is_verbose if TRUE, gbm2 will print out progress and performance of the fit.
#' 
#' @return a \code{GBMFit} object.
#' 
#' @export gbm2.fit
#' 

gbm2.fit <- function(x, y, distribution=gbm_dist("Gaussian", ...), weights=rep(1, nrow(x)), offset=rep(0, nrow(x)),
                     train_params=training_params(num_trees=100, interaction_depth=1, min_num_obs_in_node=10, 
                     shrinkage=0.001, bag_fraction=0.5, id=seq_len(nrow(x)), num_train=round(0.5 * nrow(x)), num_features=ncol(x)), 
                     var_monotone=NULL, var_names=NULL, keep_gbm_data=FALSE, cv_folds=1,
                     cv_class_stratify=FALSE, fold_id=NULL, is_verbose=FALSE) {
  
  # Check num_features makes sense
  if(train_params$num_features > ncol(x)) {
    warning("Number of features exceeds the number of predictor variables - 
            setting number of features to number of predictors")
    train_params$num_features <- ncol(x)
  }
  
  # Create gbm_data_obj
  gbm_data_obj <- gbm_data(x, y, weights, offset)
  
  # Set up groups
  distribution <- determine_groups(colnames(x), gbm_data_obj, distribution)
  
  # Set-up variable containers
  variables <- var_container(gbm_data_obj, var_monotone, var_names)
  
  # Create strata
  distribution <- create_strata(gbm_data_obj, train_params, distribution)
  
  # Process data obj and validate
  gbm_data_obj <- convert_factors(gbm_data_obj)
  gbm_data_obj <- validate_gbm_data(gbm_data_obj, distribution)
  
  # Order the data
  gbm_data_obj <- order_data(gbm_data_obj, distribution, train_params)
  
  # Get CV groups
  cv_groups <- create_cv_groups(gbm_data_obj, distribution, train_params, cv_folds,
                                cv_class_stratify, fold_id)
  # Create fitted object
  gbm_fit_obj <- gbm_cross_val(gbm_data_obj, distribution, train_params, variables,
                               cv_folds, cv_groups, is_verbose)
  
  # Keep original data
  if(keep_gbm_data) {
    gbm_fit_obj$gbm_data_obj <- gbm_data_obj
  }
  # Reorder if necessary
  gbm_fit_obj <- reorder_fit(gbm_fit_obj, distribution)
  gbm_fit_obj$fold_ids <- fold_id
  
  return(gbm_fit_obj)
} 