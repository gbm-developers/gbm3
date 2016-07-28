#' Create CV Groups
#' 
#' Methods that create cv groups for gbm2 fitting.
#' Used internally, within \code{\link{gbm_cross_val}}, to split data and distribution specific parameters into
#' train and test folds.
#' 
#' @usage create_cv_groups(gbm_data_obj, gbm_dist_obj, train_params, cv_folds, cv_class_stratify, fold_id)
#' 
#' @param gbm_data_obj a GBMData object containing all of the data used to fit a gbm model. 
#' 
#' @param gbm_dist_obj a GBMDist object specifying the distribution and any additional parameters needed.
#' 
#' @param train_params a GBMTrainParams object containing generic parameters defining how the model should be
#' trained.
#' 
#' @param cv_folds a positive integer specifying the number of folds to be used in cross-validation of the gbm
#' fit. If cv_folds > 1 then cross-validation is performed; the default of cv_folds is 1.
#' 
#' @param cv_class_stratify a bool specifying whether or not to stratify via response outcome. Currently only 
#' applies to "Bernoulli" distribution and defaults to false.
#' 
#' @param fold_id An optional vector of values identifying what fold
#' each observation is in.
#' 
#' @return a vector of numbers mapping each row of data to a cv fold.
#'  
#' @export create_cv_groups

create_cv_groups <- function(gbm_data_obj, gbm_dist_obj, train_params, cv_folds,
                             cv_class_stratify, fold_id) {
  # Initial checks
  check_if_gbm_data(gbm_data_obj)
  check_if_gbm_dist(gbm_dist_obj)
  check_if_gbm_train_params(train_params)
  check_if_natural_number(cv_folds)
  if(!is.logical(cv_class_stratify) || (length(cv_class_stratify) > 1)) {
    stop("cv_class_stratify must be a logical")
  }
  
  UseMethod("create_cv_groups", gbm_dist_obj)
}

create_cv_groups.default<- function(gbm_data_obj, gbm_dist_obj, train_params, cv_folds,
                                             cv_class_stratify, fold_id) {
  if(!is.null(fold_id)) {
    return(fold_id)
  } else {
    return( sample(rep(seq_len(cv_folds), length=train_params$num_train)) )
  }
}

create_cv_groups.BernoulliGBMDist <- function(gbm_data_obj, gbm_dist_obj, train_params, cv_folds,
                                             cv_class_stratify, fold_id) {
  if(cv_class_stratify) {
    # Number in each class
    Ones <- tabulate(gbm_data_obj$y[seq_len(train_params$num_train)])
    Zeros <- length(gbm_data_obj$y[seq_len(train_params$num_train)])-Ones
    
    smallGroup <- min(c(Ones, Zeros))
    
    if(smallGroup < cv_folds){
      stop(
        paste("The smallest class has only"
              ,smallGroup
              ,"objects in the training set. Can't do"
              ,cv_folds, "fold cross-validation.")
      )
    }
    
    cv_group <- vector(length=train_params$num_train)
    cv_group[gbm_data_obj$y[seq_len(train_params$num_train)] == 0] <- sample(rep(seq_len(cv_folds), length=Zeros))
    cv_group[gbm_data_obj$y[seq_len(train_params$num_train)] == 1] <- sample(rep(seq_len(cv_folds), length=Ones))
    return(cv_group)
    
  } else if(!is.null(fold_id)) {
    return(fold_id)
  } else {
    return( sample(rep(seq_len(cv_folds), length=train_params$num_train)) )
  }
}

create_cv_groups.PairwiseGBMDist <- function(gbm_data_obj, gbm_dist_obj, train_params, cv_folds,
                                             cv_class_stratify, fold_id) {
  # Split into CV folds at group boundaries
  samp <- sample(rep(seq_len(cv_folds), length=nlevels(gbm_dist_obj$group)))
  return(samp[as.integer(gbm_dist_obj$group[seq_len(train_params$num_train)])])
}


