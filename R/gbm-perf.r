#' GBMT Performance
#' 
#' Estimates optimal number of boosting iterations given a
#' \code{GBMFit} object and optionally plots various performance
#' measures.
#'
#' @param object a \code{GBMFit} object created from an initial
#' call to \code{\link{gbmt}} or \code{\link{gbm}}.
#'
#' @param plot.it an indicator of whether or not to plot the
#' performance measures. Setting \code{plot.it=TRUE} creates two
#' plots. The first plot plots the train error (in black)
#' and the validation error (in red) versus the iteration
#' number. The scale of the error measurement, shown on the left
#' vertical axis, depends on the \code{distribution} argument used in
#' the initial call.
#'
#' @param oobag.curve indicates whether to plot the out-of-bag
#' performance measures in a second plot.
#'
#' @param overlay if TRUE and oobag.curve=TRUE then a right y-axis is
#' added to the training and test error plot and the estimated
#' cumulative improvement in the loss function is plotted versus the
#' iteration number.
#' 
#' @param method indicate the method used to estimate the optimal
#' number of boosting iterations. \code{method="OOB"} computes the
#' out-of-bag estimate and \code{method="test"} uses the test (or
#' validation) dataset to compute an out-of-sample
#' estimate. \code{method="cv"} extracts the optimal number of
#' iterations using cross-validation if \code{gbm} was called with
#' \code{cv.folds}>1.
#' 
#' @param main the main title for the plot. Defaults to \code{main =
#' ""}.
#'
#' @return \code{gbm.perf} returns the estimated optimal number of iterations.
#' The method of computation depends on the \code{method} argument.
#' 
#' @seealso \code{\link{gbmt}} \code{\link{gbmt_performance}}
#' \code{\link{plot.GBMTPerformance}}
#'
#' @keywords nonlinear survival nonparametric tree
#'
#' @export
gbm.perf <- function(object,
                     plot.it=TRUE,
                     oobag.curve=FALSE,
                     overlay=TRUE,
                     method,
                     main="") {
    if(!is.logical(plot.it) || (length(plot.it)) > 1 || is.na(plot.it))
        stop("plot.it must be a logical - excluding NA")

    performance <- gbmt_performance(object, method)
    if (plot.it) {
        plot(performance,
             out_of_bag_curve=oobag.curve,
             overlay=overlay,
             main=main)
    }

    as.numeric(performance)
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
##'
##' @return a GBMTPerformance object, which is a number - the optimal
##' iteration number - with various attributes.
##' @export
gbmt_performance <- function(gbm_fit_obj, method) {
    check_if_gbm_fit(gbm_fit_obj)

    ## guess the method

    if (missing(method)) {
        method <- guess_error_method(gbm_fit_obj)
        message("Using ", method, " method...")
    }
    
    result <-
        switch(method,
               OOB=best_iter_out_of_bag(gbm_fit_obj),
               cv=best_iter_cv(gbm_fit_obj),
               test=best_iter_test(gbm_fit_obj),
               stop("method must be cv, test, or OOB"))

    attr(result, 'decoration') <-
        list(method=method,
             gbm_fit_obj=gbm_fit_obj)
    class(result) <- "GBMTPerformance"
    result
}


##' @export
as.double.GBMTPerformance <- function(x, ...) {
    as.double(unclass(x))
}

##' @export
print.GBMTPerformance <- function(x, ...) {
    decoration <- attr(x, 'decoration')
    method_descriptor <-
        switch(decoration$method,
               cv="cross-validation",
               test="test-set",
               OOB="out-of-bag",
               stop("Unknown method."))
    
    cat("The best ", method_descriptor, " iteration was ", x, ".\n",
        sep="")
    invisible(x)
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
    decoration <- attr(x, 'decoration')
    perf_plot(decoration$gbm_fit_obj, x,
              out_of_bag_curve, overlay,
              decoration$method,
              main)
}

#### Helper functions ####
best_iter_test <- function(gbm_fit_obj) {
  check_if_gbm_fit(gbm_fit_obj)
  best_iter_test <- which.min(iteration_error(gbm_fit_obj, 'valid'))
  return(best_iter_test)
}

best_iter_cv <- function(gbm_fit_obj) {
  check_if_gbm_fit(gbm_fit_obj)

  if(!has_cross_validation(gbm_fit_obj)) {
      stop('In order to use method="cv" gbm must be called with cv_folds>1.')
  }
  
  best_iter_cv <- which.min(iteration_error(gbm_fit_obj, 'cv'))
  return(best_iter_cv)
}

best_iter_out_of_bag <- function(gbm_fit_obj) {
  check_if_gbm_fit(gbm_fit_obj)
  if(gbm_fit_obj$params$bag_fraction==1)
    stop("Cannot compute OOB estimate or the OOB curve when bag_fraction=1")
  if(all(!is.finite(gbm_fit_obj$oobag.improve)))
    stop("Cannot compute OOB estimate or the OOB curve. No finite OOB estimates of improvement")
  
  message("OOB generally underestimates the optimal number of iterations although predictive performance is reasonably competitive.
            Using cv_folds>1 when calling gbm usually results in improved predictive performance.")
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
