# Num Train Obs
# Function that generates the number of training
# observations - this function is introduced for 
# backwards compatibility with train.fraction

num_train_obs <- function(train.fraction, gbm_train_params, distribution_obj) {
  check_if_gbm_dist(distribution_obj)
  check_if_gbm_params(gbm_train_params)
  
  if((length(train.fraction) > 1) || !is.double(train.fraction)
     || is.null(train.fraction) || (train.fraction > 1.0)
     || (train.fraction < 0.0)) {
     stop("train.fraction must be a double between 0 and 1.")
  }
  if(distribution_obj$name != "Pairwise") {
    gbm_train_params$num_train <- floor(train.fraction * length(unique(gbm_train_params$id)))
  } else {
    # Split into train and validation set, at group boundary
    num.groups.train <- max(1, round(train.fraction * nlevels(distribution_obj$group)))
    
    # include all groups up to the num.groups.train
    gbm_train_params$num_train <- max(which(distribution_obj$group==levels(distribution_obj$group)[num.groups.train]))
  }
  
  return(gbm_train_params)
}