#' Extract observations in fold
#' 
#' Extract the relevant observations to fit a CV fold - reorders appropriately.
#' 
#' @usage extract_obs_in_fold(gbm_data_obj, gbm_dist_obj, train_params, cv_groups, fold_num, is_verbose)
#' 
#' @param gbm_data_obj a GBMData object containing all of the data used to fit a gbm model. 
#' 
#' @param gbm_dist_obj a GBMDist object specifying the distribution and any additional parameters needed.
#' 
#' @param train_params a GBMTrainParams object containing generic parameters defining how the model should be
#' trained.
#' 
#' @param cv_groups vector of integers specifying which folds each row of data belongs to.
#' 
#' @param fold_num integer (>=1) specifying which fold under consideration.
#'
#' @param is_verbose if TRUE, will print out progress and performance of the fitting.
#' 
#' @return list of input gbm objects updated so cross_validated data is ready for gbm_call
#' 
#' @export

extract_obs_in_fold <- function(gbm_data_obj, gbm_dist_obj, train_params, cv_groups, fold_num) {
  # Observations in the training set
  obs_in_training_set <- train_params$id %in% as.integer(names(train_params$num_rows_per_obs[seq_len(train_params$num_train)]))
  train_params$id <- train_params$id[obs_in_training_set]
  
  # Observations in cv_group
  obs_id_in_cv_group <- train_params$id %in% as.integer(names(train_params$num_rows_per_obs[seq_len(train_params$num_train)]))[(cv_groups == fold_num)]

  # Extract relevent data - split into training and validation sets 
  # Calculate new number of training rows and observations
  train_params$num_train <- length(unique(train_params$id[!obs_id_in_cv_group]))
  train_params$num_train_rows <- length(train_params$id[!obs_id_in_cv_group])
  #train_params$id <- train_params$id[!obs_id_in_cv_group]
  
  gbm_data_obj <- split_and_join(gbm_data_obj, train_params, obs_in_training_set, obs_id_in_cv_group)
  gbm_dist_obj <- update_fold_dist_data(gbm_dist_obj, gbm_data_obj, train_params, obs_in_training_set, obs_id_in_cv_group)
  
  # NB: Update train_params id afterwards
  train_params$id <- train_params$id[!obs_id_in_cv_group]
  
  return(list("data"=gbm_data_obj, "dist"=gbm_dist_obj, "params"=train_params))
}
