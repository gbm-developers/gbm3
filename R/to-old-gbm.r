#' Convert GBMFit to previous gbm object
#'
#' Function that takes a \code{GBMFit} object produced from a call
#' to \code{\link{gbmt}} and converts it to the form of a gbm object from
#' Version 2.1 of the package.
#'
#' @usage  to_old_gbm(gbm_fit_obj)
#'
#' @param gbm_fit_obj a \code{GBMFit} object produced by a call to \code{\link{gbmt}}.
#'
#' @return a \code{gbm} object of the form from Version 2.1 of the package.
#'
#' @author James Hickey
#'
#' @export

to_old_gbm <- function(gbm_fit_obj) {
  # Check input
  check_if_gbm_fit(gbm_fit_obj)

  # Convert gbm_fit_obj to old API
  gbm_fit_old <- list()
  gbm_fit_old$initF <- gbm_fit_obj$initF
  gbm_fit_old$train.error <- iteration_error(gbm_fit_obj, 'train')
  gbm_fit_old$valid.error <- iteration_error(gbm_fit_obj, 'valid')
  gbm_fit_old$trees <- trees(gbm_fit_obj)
  gbm_fit_old$c.splits <- gbm_fit_obj$c.splits
  gbm_fit_old$oobag.improve <- gbm_fit_obj$oobag.improve
  gbm_fit_old$fit <- gbm_fit_obj$fit

  gbm_fit_old$bag.fraction <- gbm_fit_obj$params$bag_fraction
  gbm_fit_old$distribution <- list(name=tolower(distribution_name(gbm_fit_obj)))
  gbm_fit_old$interaction.depth <- gbm_fit_obj$params$interaction_depth
  gbm_fit_old$n.minobsinnode <- gbm_fit_obj$params$min_num_obs_in_node
  gbm_fit_old$n.trees <- length(gbm_fit_old$trees)
  gbm_fit_old$nTrain <- gbm_fit_obj$params$num_train_rows
  gbm_fit_old$nTrainPats <- gbm_fit_obj$params$num_train
  gbm_fit_old$patient.id <- gbm_fit_obj$params$id
  gbm_fit_old$mFeatures <- gbm_fit_obj$params$num_features
  gbm_fit_old$train.fraction <- gbm_fit_obj$params$train_fraction
  gbm_fit_old$response.name <- gbm_fit_obj$response_name
  gbm_fit_old$shrinkage <- gbm_fit_obj$params$shrinkage
  gbm_fit_old$var.levels <- gbm_fit_obj$variables$var_levels
  gbm_fit_old$var.monotone <- gbm_fit_obj$variables$var_monotone
  gbm_fit_old$var.names <- gbm_fit_obj$variables$var_names
  gbm_fit_old$var.type <- gbm_fit_obj$variables$var_type
  gbm_fit_old$verbose <- gbm_fit_obj$is_verbose
  gbm_fit_old$strata <- gbm_fit_obj$distribution$strata
  gbm_fit_old$sorted <- gbm_fit_obj$distribution$sorted
  gbm_fit_old$prior.node.coeff.var <- ifelse(is.null(gbm_fit_obj$distribution$prior_node_coeff_var), 1000,
                                         gbm_fit_obj$distribution$prior_node_coeff_var)

  if(!is.null(gbm_fit_obj$gbm_data_obj)) {
    # put the observations back in - these are ordered according to id and group
    data <- gbm_fit_obj$gbm_data_obj
    if(distribution_name(gbm_fit_obj) == "CoxPH") {
      gbm_fit_old$data <- list(y=data$y, x=data$x,x.order=data$x_order, offset=data$offset,
                               Misc=unlist(get_misc(gbm_fit_obj$distribution)), w=data$weights,
                               i.order=gbm_fit_obj$distribution$time_order)
    } else {
      gbm_fit_old$data <- list(y=data$y, x=data$x,x.order=data$x_order, offset=data$offset,
                               Misc=unlist(get_misc(gbm_fit_obj$distribution)), w=data$weights)
    }

  } else {
    gbm_fit_old$data <- NULL
  }
  gbm_fit_old$tied.times.methods <- gbm_fit_obj$distribution$ties
  gbm_fit_old$ord.group <- gbm_fit_obj$distribution$group_order
  gbm_fit_old$cv.folds <- gbm_fit_obj$cv_folds
  gbm_fit_old$cv.error <- iteration_error(gbm_fit_obj, 'cv')
  gbm_fit_old$cv.fitted <- gbm_fit_obj$cv_fitted
  gbm_fit_old$Terms <- gbm_fit_obj$Terms
  gbm_fit_old$call <- gbm_fit_obj$call
  gbm_fit_old$m <- gbm_fit_obj$m
  gbm_fit_old$num.classes <- gbm_fit_obj$num.classes


  class(gbm_fit_old) <- "gbm"

  message("Converted to old gbm object - this will not work with new packages functions")

  return(gbm_fit_old)

}
