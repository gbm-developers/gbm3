#' GBM2 Fit
#' 
#' Wrapper function for calling the C++ gbm function
#' 
#' @usage gbm_fit(gbm_data_obj, gbm_dist_obj, train_params, var_container, is_verbose)
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
#' @param is_verbose if TRUE, will print out progress and performance of the fitting.
#'
#' @return a fitted gbm object

gbm_fit <- function(gbm_data_obj, gbm_dist_obj, train_params, var_container, is_verbose) {
  # Check inputs
  check_if_gbm_data(gbm_data_obj)
  check_if_gbm_dist(gbm_dist_obj)
  check_if_gbm_train_params(train_params)
  check_if_gbm_var_container(var_container)
  
  fit <- .Call("gbm",
                Y=as.matrix(as.data.frame(gbm_data_obj$y)),
                Offset=as.double(gbm_data_obj$offset),
                X=as.matrix(as.data.frame(gbm_data_obj$x)),
                X.order=as.integer(gbm_data_obj$x_order),
                sorted=as.matrix(as.data.frame(gbm_dist_obj$sorted)),
                Strata = as.integer(gbm_dist_obj$strata),
                weights=as.double(gbm_data_obj$weights),
                Misc=get_misc(gbm_dist_obj),
                prior.node.coeff.var = ifelse(is.null(gbm_dist_obj$prior_node_coeff_var), as.double(0),
                                                         as.double(gbm_dist_obj$prior_node_coeff_var)),
                id = as.integer(train_params$id),
                var.type=as.integer(var_container$var_type),
                var.monotone=as.integer(var_container$var_monotone),
                distribution=as.character(tolower(gbm_dist_obj$name)),
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
  return(fit)
}