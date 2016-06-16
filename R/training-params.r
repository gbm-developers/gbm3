#' Training parameters
#' 
#' Class that contains the training parameters for the gbm model
#' 
#' @usage training_params(num_trees=100, interaction_depth=1,
#'  min_num_obs_in_node=10, shrinkage=0.001, bag_fraction=0.5, id,  num_train)
#' 
#' @param num_trees Number of trees used in the fit.
#' 
#' @param interaction_depth Maximum depth of each tree
#' 
#' @param min_num_obs_in_node Minimum number of observations each
#' node in a tree must have.
#' 
#' @param shrinkage shrinkage parameter applied to each tree in the expansion. 
#' Also known as the learning rate or step-size reduction.
#' 
#' @param bag_fraction fraction of independent training observations selected to
#' create the next tree in the expansion.  Introduces randomness in the model fit; 
#' if bag_fraction < 1 then running the same model twice  will result in similar but different
#' fits.
#' 
#' @param id optional vector of integers, specifying which rows in the data correspond
#' to which observations. Individual observations may have many rows of data associated with them.
#' 
#' @param num_train number of unique observations used in training the model.
#' 
#' @param num_features number of random features/columns to use in training model
#' 
#' @return training parameters object
#' 
#' @export

training_params <- function(num_trees=100, interaction_depth=1,
                            min_num_obs_in_node=10, shrinkage=0.001, bag_fraction=0.5,
                            id, num_train, num_features) {
  # Check the parameters
  check_if_natural_number(num_trees, "number of trees")
  check_if_natural_number(interaction_depth, "interaction depth")
  check_if_natural_number(min_num_obs_in_node, "minimum number of node obs")
  check_if_natural_number(num_train, "number of training observations")
  check_if_natural_number(num_features, "number of features")
  checkID(interaction_depth)
  
  if(length(id) < num_train) {
    stop("Number of unique observation ids is less than the amount of training data")
  } else if(!is.atomic(id) || !all(id == as.integer(id))) {
    stop("Observation ids must be a vector of inter")
  }
  
  if(!is.double(shrinkage) || (shrinkage > 1.0) 
     || is.infinite(shrinkage) || (length(shrinkage) > 1)) {
    stop("Shrinkage should be a double < 1.0")
  }
  if(!is.double(bag_fraction) || (bag_fraction > 1.0)
     || (length(bag_fraction) > 1) || (bag_fraction < 0.0)) {
    stop("Bag fraction should be a double between 0.0 and 1.0")
  }
  
  object <- structure(list("num_trees"=num_trees, "interaction_depth"=interaction_depth,
                 "min_num_obs_in_node"=min_num_obs_in_node, "shrinkage"=shrinkage,
                 "bag_fraction"=bag_fraction, "id"= id, "num_train"=num_train,
                 "num_features"=num_features),
                 class="trainingParams")
  
  return(object)

}

