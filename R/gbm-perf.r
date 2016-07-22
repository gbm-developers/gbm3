#' GBM2 Performace
#' 
#' Estimates optimal number of boosting iterations given a \code{GBMFit} object and
#' optionally plots various performance measures.
#' 
#' @usage gbm_perf(gbm_fit_obj, plot_it=TRUE, out_of_bag_curve=FALSE, overlay=TRUE, method, main="")
#' 
#' @param gbm_fit_obj a \code{GBMFit} created from an initial call to
#' \code{\link{gbm2}}.
#' 
#' @param plot_it an indicator of whether or not to plot the performance
#' measures. Setting \code{plot_it=TRUE} creates two plots. The first plot
#' plots \code{gbm_fit_obj$train.error} (in black) and \code{gbm_fit_obj$valid.error} (in
#' red) versus the iteration number. The scale of the error measurement, shown
#' on the left vertical axis, depends on the \code{distribution} argument used
#' in the initial call to \code{\link{gbm2}}.
#' 
#' @param out_of_bag_curve indicates whether to plot the out-of-bag performance
#' measures in a second plot.
#' 
#' @param overlay if TRUE and out_of_bag_curve=TRUE then a right y-axis is added to
#' the training and test error plot and the estimated cumulative improvement in
#' the loss function is plotted versus the iteration number.
#' 
#' @param method indicate the method used to estimate the optimal number of
#' boosting iterations. \code{method="OOB"} computes the out-of-bag estimate
#' and \code{method="test"} uses the test (or validation) dataset to compute an
#' out-of-sample estimate. \code{method="cv"} extracts the optimal number of
#' iterations using cross-validation if \code{gbm2} was called with \code{cv_folds}>1.
#' 
#' @param main the main title for the plot. Defaults to \code{main = ""}.
#' 
#' @return \code{gbm_perf} returns the estimated optimal number of iterations.
#' The method of computation depends on the \code{method} argument.
#' @seealso \code{\link{gbm2}}
#' @keywords nonlinear survival nonparametric tree
#' @export 
#' 


gbm_perf <- function(gbm_fit_obj, plot_it=TRUE, 
                                        out_of_bag_curve=FALSE,
                                        overlay=TRUE,
                                        method,
                                        main="") {
  # Initial checks
  check_if_gbm_fit(gbm_fit_obj)
  if ( missing( method ) )
    stop("requires method parameter to determine performance")
  
  if(!is.element(method,c("OOB","test","cv")))
    stop("method must be cv, test, or OOB")
  
  if(!is.logical(plot_it) || (length(plot_it)) > 1)
    stop("plot_it must be a logical")
  
  if(!is.logical(plot_it) || (length(plot_it)) > 1)
    stop("plot_it must be a logical")
  
  best_iter <- switch(method,
                      OOB=best_iter_out_of_bag(gbm_fit_obj),
                      cv=best_iter_cv(gbm_fit_obj),
                      test=best_iter_test(gbm_fit_obj))
  if(plot_it)
    perf_plot(gbm_fit_obj, best_iter, out_of_bag_curve, overlay, method, main)
  
  return(best_iter)  
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
    stop("In order to use method=\"cv\" gbm must be called with cv.folds>1.")
    warning("cross-validation error is not computed for any additional iterations run using gbm.more().")
  best_iter_cv <- which.min(gbm_fit_obj$cv_error)
  return(best_iter_cv)
}

best_iter_out_of_bag <- function(gbm_fit_obj) {
  check_if_gbm_fit(gbm_fit_obj)
  if(gbm_fit_obj$params$bag_fraction==1)
    stop("Cannot compute OOB estimate or the OOB curve when bag.fraction=1")
  if(all(!is.finite(gbm_fit_obj$oobag.improve)))
    stop("Cannot compute OOB estimate or the OOB curve. No finite OOB estimates of improvement")
  
  warning("OOB generally underestimates the optimal number of iterations although predictive performance is reasonably competitive.
            Using cv.folds>0 when calling gbm usually results in improved predictive performance.")
  smoother <- generate_smoother_oobag(gbm_fit_obj)
  best_iter_oob <- x[which.min(-cumsum(smoother$y))]
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