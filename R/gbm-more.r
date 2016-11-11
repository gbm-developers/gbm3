#' Perform additional boosting
#' 
#' Method to perform additional boosting using a \code{GBMFit} object
#' - does not support further cross validation.
#' 
#' @param gbm_fit_obj a \code{GBMFit} object produced using
#' \code{\link{gbmt}}.  This object describes the boosted model on
#' which to perform additional boosting.
#' 
#' @param num_new_trees a positive integer specifying how many
#' additional iterations to perform.  This has a default value of
#' \code{100}.
#' 
#' @param data a \code{data.frame} or \code{matrix} containing the new
#' values for the predictor and response variables for the additional
#' iterations.  The names of the variables must match those appearing
#' in the original fit (as well as the number of rows), and this value
#' defaults to \code{NULL}.  With a value of \code{NULL} the original
#' data may be used for the additional boosting, if no original or new
#' data is specified an error will be thrown.
#' 
#' @param weights an atomic vector of doubles specifying the
#' importance of each row of the \code{data} in the additional
#' iterations. If the previous data used is kept within
#' \code{gbm_fit_obj}; then the weights are extracted from the stored
#' \code{GBMData} object.
#' 
#' @param offset an atomic vector of doubles specifying the offset for
#' each response value in the data used for additional boosting.
#' 
#' @param is_verbose a logical specifying whether or not the
#' additional fitting should run "noisely" with feedback on progress
#' provided to the user.
#' 
#' @return the input \code{GBMFit} object with additional iterations
#' provided for the fit.
#' 
#' @export 
gbm_more <- function(gbm_fit_obj, num_new_trees=100,
                     data=NULL, weights=NULL, offset=NULL,
                     is_verbose=FALSE){
  the_call <- match.call()
  
  # Check inputs
  check_if_gbm_fit(gbm_fit_obj)
  check_if_natural_number(num_new_trees)

  if(!is.logical(is_verbose) || (length(is_verbose) > 1) || is.na(is_verbose)) {
      stop("is_verbose must be a logical - not NA")
  }
  
  if(is.null(gbm_fit_obj$gbm_data_obj) && is.null(data)) {
      stop("keep_data was set to FALSE on original gbmt call and argument 'data' is NULL")
  }

  if (has_cross_validation(gbm_fit_obj)) {
      warning("gbm_more is incompatible with cross-validation; losing cv results.")
  }
  
  # If no data is stored create appropriate data objects 
  if(is.null(gbm_fit_obj$gbm_data_obj)) {
    m <- eval(gbm_fit_obj$model, parent.frame())
    Terms <- attr(m, "terms")
    a <- attributes(Terms)
    y <- as.vector(model.extract(m, "response"))
    offset <- model.extract(m,offset)
    x <- model.frame(delete.response(Terms),
                     data,
                     na.action=na.pass)
    w <- weights
    
    # If offset or weights NULL change to simple values
    if(length(weights)==0) w <- rep(1, nrow(x))
    if(length(offset)==0) offset <- rep(0, nrow(x))
    
    # Create gbm_data
    gbm_data_obj <- gbm_data(x, y, weights, offset)
    
    # Get distribution and set up groups
    distribution <- determine_groups(colnames(data), gbm_data_obj, gbm_fit_obj$distribution)
    
    # Create strata
    distribution <- create_strata(gbm_data_obj, gbm_fit_obj$params, distribution)
    
    # Process data obj and validate
    gbm_data_obj <- convert_factors(gbm_data_obj)
    gbm_data_obj <- validate_gbm_data(gbm_data_obj, distribution)
    
    # Order the data
    gbm_data_obj <- order_data(gbm_data_obj, distribution, gbm_fit_obj$params)
    
  } else {
    gbm_data_obj <- gbm_fit_obj$gbm_data_obj
    distribution <- gbm_fit_obj$distribution
  }
  
  # Reorder fit so as to be in order of data
  if (distribution_name(gbm_fit_obj) == "Pairwise") 
  {
    gbm_fit_obj$fit   <- gbm_fit_obj$fit[gbm_fit_obj$distribution$group_order] # object$fit is stored in the original order
  } else if (distribution_name(gbm_fit_obj) == "CoxPH") {
    gbm_fit_obj$fit <- gbm_fit_obj$fit[gbm_fit_obj$time_order]
  }

  # Call GBM package
  gbm_more_fit <-
      .Call("gbm",
            Y=as.matrix(as.data.frame(gbm_data_obj$y)),
            intResponse = as.matrix(cbind(distribution$strata, distribution$sorted)),
            Offset=as.double(gbm_data_obj$offset),
            X=as.matrix(as.data.frame(gbm_data_obj$x)),
            X.order=as.integer(gbm_data_obj$x_order),
            weights=as.double(gbm_data_obj$weights),
            Misc=get_misc(distribution),
            prior.node.coeff.var = ifelse(is.null(distribution$prior_node_coeff_var), as.double(0),
                as.double(distribution$prior_node_coeff_var)),
            id = as.integer(gbm_fit_obj$params$id),
            var.type=as.integer(gbm_fit_obj$variables$var_type),
            var.monotone=as.integer(gbm_fit_obj$variables$var_monotone),
            distribution=gbm_call_dist_name(distribution),
            n.trees=as.integer(num_new_trees),
            interaction.depth=as.integer(gbm_fit_obj$params$interaction_depth),
            n.minobsinnode=as.integer(gbm_fit_obj$params$min_num_obs_in_node),
            shrinkage=as.double(gbm_fit_obj$params$shrinkage),
            bag.fraction=as.double(gbm_fit_obj$params$bag_fraction),
            nTrainRows=as.integer(gbm_fit_obj$params$num_train_rows),
            nTrainObs = as.integer(gbm_fit_obj$params$num_train),
            mFeatures=as.integer(gbm_fit_obj$params$num_features),
            fit.old=as.double(gbm_fit_obj$fit),
            n.cat.splits.old=as.integer(length(gbm_fit_obj$c.splits)),
            n.trees.old=length(trees(gbm_fit_obj)),
            par_details=gbm_fit_obj$par_details,
            verbose=as.integer(is_verbose),
            PACKAGE = "gbm")
  
  # Set class
  class(gbm_more_fit) <- "GBMFit"
  
  # Store correct parameters
  gbm_more_fit$distribution <- distribution
  gbm_more_fit$params <- gbm_fit_obj$params
  gbm_more_fit$variables <- gbm_fit_obj$variables
  gbm_more_fit$Terms     <- gbm_fit_obj$Terms
  
  # Transfer old results across
  gbm_more_fit$initF         <- gbm_fit_obj$initF
  gbm_more_fit$train.error   <- c(gbm_fit_obj$train.error,
                                  gbm_more_fit$train.error)
  gbm_more_fit$valid.error   <- c(gbm_fit_obj$valid.error,
                                  gbm_more_fit$valid.error)
  gbm_more_fit$oobag.improve <- c(gbm_fit_obj$oobag.improve,
                                  gbm_more_fit$oobag.improve)
  gbm_more_fit$trees         <- c(trees(gbm_fit_obj),
                                  gbm_more_fit$trees)
  gbm_more_fit$c.splits      <- c(gbm_fit_obj$c.splits,
                                  gbm_more_fit$c.splits)
  gbm_more_fit$params$num_trees <- length(gbm_more_fit$trees)
    
  # Compatibility with 1.6
  gbm_more_fit$num.classes <- gbm_fit_obj$num.classes

  # Reorder if necessary 
  gbm_more_fit <- reorder_fit(gbm_more_fit, distribution)

  # Store data 
  if(!is.null(gbm_fit_obj$gbm_data_obj)) {
    gbm_more_fit$gbm_data_obj <- gbm_fit_obj$gbm_data_obj
  }
  gbm_more_fit$model <- gbm_fit_obj$model
  gbm_more_fit$par_details <- gbm_fit_obj$par_details
  gbm_more_fit$call <- the_call
  
  return(gbm_more_fit)
  
}
