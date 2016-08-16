#' GBMT Performance
#' 
#' Estimates optimal number of boosting iterations given a
#' \code{GBMFit} object and optionally plots various performance
#' measures.
#'
#' @param plot_it an indicator of whether or not to plot the
#' performance measures. 
#' @inheritParams gbmt_performance
#' @inheritParams plot.GBMTPerformance
#'
#' @return \code{gbm_perf} returns the estimated optimal number of
#' iterations.  The method of computation depends on the \code{method}
#' argument.
#' @seealso \code{\link{gbmt}} \code{\link{gbmt_performance}}
#' \code{\link{plot.GBMTPerformance}}
#' @keywords nonlinear survival nonparametric tree
#' @export
gbm_perf <- function(gbm_fit_obj, plot_it=TRUE, 
                     out_of_bag_curve=FALSE,
                     overlay=TRUE,
                     method,
                     main="") {
    if(!is.logical(plot_it) || (length(plot_it)) > 1 || is.na(plot_it))
        stop("plot_it must be a logical - excluding NA")

    ## guess the method (to match old gbm.perf)
    if (missing(method)) {
        method <- guess_error_method(gbm_fit_obj)
        message("Using ", method, " method...")
    }

    performance <- gbmt_performance(gbm_fit_obj, method)
    if (plot_it) {
        plot(performance,
             out_of_bag_curve=out_of_bag_curve,
             overlay=overlay,
             main=main)
    }

    summary(performance)
}

##' Get performance details for gbm fit
##'
##' \code{gbmt_performance} estimates the optimal number of boosting
##' iterations from a model fit by \code{\link{gbmt}}.  The precise
##' method used depends on the \code{method} parameter.
##' 
##' @param gbm_fit_obj a \code{GBMFit} created from an initial call to
##' \code{\link{gbmt}}.
##'
##' @param method indicate the method used to estimate the optimal
##' number of boosting iterations. \code{method="OOB"} computes the
##' out-of-bag estimate and \code{method="test"} uses the test (or
##' validation) dataset to compute an out-of-sample
##' estimate. \code{method="cv"} extracts the optimal number of
##' iterations using cross-validation if \code{gbmt} was called with
##' \code{cv_folds}>1.
##' @return a GBMTPerformance object
##' @export
gbmt_performance <- function(gbm_fit_obj, method) {
    check_if_gbm_fit(gbm_fit_obj)
    
    best_iter <-
        switch(method,
               OOB=best_iter_out_of_bag(gbm_fit_obj),
               cv=best_iter_cv(gbm_fit_obj),
               test=best_iter_test(gbm_fit_obj),
               stop("method must be cv, test, or OOB"))

    result <- list(best_iter=best_iter,
                   method=method,
                   gbm_fit_obj=gbm_fit_obj)
    class(result) <- "GBMTPerformance"
    result
}

##' @export
summary.GBMTPerformance <- function(object, ...) {
    object$best_iter
}

##' Plot GBM performance details
##'
##' The train and validation error (in black and red respectively) are
##' plotted against the iteration number.  If the initial model was
##' built with cross-validation, the cross-validation error is shown
##' in green.
##'
##' The scale of the error measurement, shown on the
##' left vertical axis, depends on the \code{distribution} argument
##' used in the initial call to \code{\link{gbmt}}.
##'
##' @param x a \code{GBMTPerformance} object (created by
##' \code{\link{gbmt_performance}})
##' 
##' @param out_of_bag_curve indicates whether to plot the out-of-bag
##' performance measures in a second plot.
##' 
##' @param overlay if TRUE and out_of_bag_curve=TRUE then a right
##' y-axis is added to the training and test error plot and the
##' estimated cumulative improvement in the loss function is plotted
##' versus the iteration number.
##'
##' @param main the main title for the plot.
##'
##' @param \dots currently ignored
##'
##' @export
plot.GBMTPerformance <- function(x,
                                 out_of_bag_curve=FALSE,
                                 overlay=TRUE,
                                 main="", ...) {
    perf_plot(x$gbm_fit_obj, x$best_iter,
              out_of_bag_curve, overlay,
              x$method,
              main)
}

#### Helper functions ####
best_iter_test <- function(gbm_fit_obj) {
  check_if_gbm_fit(gbm_fit_obj)
  best_iter_test <- which.min(gbm_fit_obj$valid.error)
  return(best_iter_test)
}

best_iter_cv <- function(gbm_fit_obj) {
  check_if_gbm_fit(gbm_fit_obj)
  if(is.null(gbm_fit_obj$cv_error))
    stop("In order to use method=\"cv\" gbm must be called with cv_folds>1.")
    message("cross-validation error is not computed for any additional iterations run using gbm_more().")
  best_iter_cv <- which.min(gbm_fit_obj$cv_error)
  return(best_iter_cv)
}

best_iter_out_of_bag <- function(gbm_fit_obj) {
  check_if_gbm_fit(gbm_fit_obj)
  if(gbm_fit_obj$params$bag_fraction==1)
    stop("Cannot compute OOB estimate or the OOB curve when bag_fraction=1")
  if(all(!is.finite(gbm_fit_obj$oobag.improve)))
    stop("Cannot compute OOB estimate or the OOB curve. No finite OOB estimates of improvement")
  
  message("OOB generally underestimates the optimal number of iterations although predictive performance is reasonably competitive.
            Using cv_folds>0 when calling gbm usually results in improved predictive performance.")
  smoother <- generate_smoother_oobag(gbm_fit_obj)
  best_iter_oob <- smoother$x[which.min(-cumsum(smoother$y))]
  return(best_iter_oob)
}

generate_smoother_oobag <- function(gbm_fit_obj) {
  check_if_gbm_fit(gbm_fit_obj)
  smoother <- NULL
  x <- seq_len(gbm_fit_obj$params$num_trees)
  smoother <- loess(gbm_fit_obj$oobag.improve~x,
                    enp.target=min(max(4,length(x)/10),50))
  smoother$y <- smoother$fitted
  smoother$x <- x
  return(smoother)
}

## What error method should we use?
guess_error_method <- function(gbm_fit_obj) {
    if (has_train_test_split(gbm_fit_obj)) {
        "test"
    } else if (has_cross_validation(gbm_fit_obj)) {
        "cv"
    } else {
        "OOB"
    }
}
