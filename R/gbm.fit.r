#' @describeIn gbm Core fitting code, for experts only.
#' @export 
gbm.fit <- function(x, y,
                    offset = NULL,
                    distribution = "Bernoulli",
                    w = NULL,
                    var.monotone = NULL,
                    n.trees = 100,
                    interaction.depth = 1,
                    n.minobsinnode = 10,
                    shrinkage = 0.001,
                    bag.fraction = 0.5,
                    nTrain = NULL,
                    train.fraction = NULL,
                    mFeatures = NULL,
                    keep.data = TRUE,
                    verbose = TRUE,
                    var.names = NULL,
                    response.name = "y",
                    group = NULL,
                    tied.times.method="efron",
                    prior.node.coeff.var = 1000,
                    strata = NA, obs.id = 1:nrow(x)) {
  # Highlight new API
  warning("gbm is depracated - using gbmt_fit...")
  
  # Reconstruct data
  if(is.null(var.names)) {
    var.names <- get_var_names(x)
  }
 
  # Calculate nTrain if necessary 
  if(!is.null(nTrain) && !is.null(train.fraction)) {
    stop("Parameters 'nTrain' and 'train.fraction' cannot both be specified")
    
  } else if(!is.null(train.fraction)) {
    warning("Parameter 'train.fraction' of gbm.fit is deprecated, please specify 'nTrain' instead")
    nTrain <- floor(train.fraction*length(unique(obs.id)))
    
  } else if(is.null(nTrain)) {
    # both undefined, use all training data
    nTrain <- length(unique(obs.id))
  }
  
  # Set offset/ weights if not specified
  if(is.null(w)) {
    weights <- rep(1, length(obs.id))
  } 
  if(is.null(offset)){
    offset <- rep(0, length(obs.id))
  }
  
  # Set distribution object - put in all possible additional parameters (this will generate warnings)
  if (is.character(distribution)){ distribution <- list(name=distribution) }
  dist_obj <- create_dist_obj_for_gbmt_fit(distribution, tied.times.method, strata, prior.node.coeff.var)
  
  # Set up training parameters
  if(is.null(mFeatures)) mFeatures <- ncol(x)
  params <- training_params(num_trees=n.trees, interaction_depth=interaction.depth, min_num_obs_in_node=n.minobsinnode, 
                            shrinkage=shrinkage, bag_fraction=bag.fraction, id=obs.id,
                            num_train=nTrain, num_features=mFeatures)
  
  # Call gbm2 - formula defined from the data
  gbm_fit_obj <- gbmt_fit(x, y, distribution=dist_obj, weights=w, offset=offset,
                      train_params= params, var_monotone=var.monotone, var_names=var.names,
                      keep_gbm_data=keep.data, is_verbose=verbose)
  
  return(gbm_fit_obj)
}
