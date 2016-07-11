#' GBM2
#' 
#' Fits generalized boosted regression models - new API.
#' 
#' @usage  gbm2(formula, distribution=gbm_dist("Gaussian", ...), data, weights, offset,
#' train_params=training_params(num_trees=100, interaction_depth=1, min_num_obs_in_node=10, 
#' shrinkage=0.001, bag_fraction=0.5, id=seq_len(nrow(data)), num_train=round(0.5 * nrow(data)), num_features) 
#' var_monotone=NULL, var_names=NULL, keep_gbm_data=FALSE, is_verbose=FALSE)
#' 
#' @param formula a symbolic description of the model to be fit.  The formula may include
#' an offset term (e.g. y~offset(n) + x).
#' 
#' @param  distribution a GBMDist object specifying the distribution and any additional parameters needed.
#' 
#' @param data a data frame containing the variables in the model.  By default, the variables are taken from the 
#' environment. 
#' 
#' @param weights optional vector of weights used in the fitting process.  These weights must be positive but 
#' need not be normalized.
#' 
#' @param offset  optional vector specifying the model offset; must be positive.
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
#' @param is_verbose if TRUE, gbm2 will print out progress and performance of the fit.
#' 
#' @return a gbm2 object.
#' 
#' @export gbm2
#' 

gbm2 <- function(formula, distribution=gbm_dist("Gaussian", ...), data, weights, offset,
                 train_params=training_params(num_trees=100, interaction_depth=1, min_num_obs_in_node=10, 
                 shrinkage=0.001, bag_fraction=0.5, id=seq(nrow(data)), num_train=1, num_features=ncol(data)-1), num_train=round(0.5 * nrow(data)), num_features, 
                 var_monotone=NULL, var_names=NULL, keep_gbm_data=FALSE, is_verbose=FALSE) {
  theCall <- match.call()
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "weights", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.pass
  mf[[1]] <- as.name("model.frame")
  m <- mf
  mf <- eval(mf, parent.frame())
  Terms <- attr(mf, "terms")
  y <- model.response(mf)
  w <- model.weights(mf)
  offset <- model.offset(mf)
  
  # get the character name of the response variable
  response_name <- as.character(formula[[2]])
  
  var_names <- attributes(Terms)$term.labels
  x <- model.frame(terms(reformulate(var_names)),
                   data,
                   na.action=na.pass)

  # Create gbm_data_obj
  gbm_data_obj <- gbm_data(x, y, weights, offset)
  
  # Set up groups
  distribution <- determine_groups(colnames(data), gbm_data_obj, distribution)
  
  # Set-up variable containers
  variables <- var_container(gbm_data_obj, var_monotone, var_names)
  
  # Process data obj and validate
  gbm_data_obj <- convert_factors(gbm_data_obj)
  gbm_data_obj <- validate_gbm_data(gbm_data_obj, distribution)

  # Create strata
  distribution <- create_strata(gbm_data_obj, train_params, distribution)
  
  # Order the data
  gbm_data_obj <- order_data(gbm_data_obj, distribution, train_params)
  
  # Call the method
  gbm_fit <- .Call("gbm",
                   Y=as.matrix(as.data.frame(gbm_data_obj$y)),
                   Offset=as.double(gbm_data_obj$offset),
                   X=as.matrix(as.data.frame(gbm_data_obj$x)),
                   X.order=as.integer(gbm_data_obj$x_order),
                   sorted=as.matrix(as.data.frame(distribution$sorted)),
                   Strata = as.integer(distribution$strata),
                   weights=as.double(gbm_data_obj$weights),
                   Misc=get_misc(distribution),
                   prior.node.coeff.var = ifelse(is.null(distribution$prior_node_coeff_var), as.double(0),
                                                 as.double(distribution$prior_node_coeff_var)),
                   id = as.integer(train_params$id),
                   var.type=as.integer(variables$var_type),
                   var.monotone=as.integer(variables$var_monotone),
                   distribution=as.character(tolower(distribution$name)),
                   n.trees=as.integer(train_params$num_trees),
                   interaction.depth=as.integer(train_params$interaction_depth),
                   n.minobsinnode=as.integer(train_params$min_num_obs_in_node),
                   shrinkage=as.double(train_params$shrinkage),
                   bag.fraction=as.double(train_params$bag_fraction),
                   nTrainRows=as.integer(train_params$num_train),
                   nTrainObs = as.integer(length(unique(train_params$id[seq_len(train_params$num_train)]))),
                   mFeatures=as.integer(train_params$num_features),
                   fit.old=as.double(NA),
                   n.cat.splits.old=as.integer(0),
                   n.trees.old=as.integer(0),
                   verbose=as.integer(is_verbose),
                   PACKAGE = "gbm")
  class(gbm_fit) <- "GBMFit"
  
  # Build up extra pieces
  gbm_fit$distribution <- distribution
  gbm_fit$params <- train_params
  gbm_fit$variables <- variables
  gbm_fit$Terms <- Terms
  
  if(keep_gbm_data) {
    gbm_fit$data <- gbm_data_obj
  }
  
  # Reorder if necessary
  gbm_fit <- reorder_fit(gbm_fit, distribution)
  return(gbm_fit)
}