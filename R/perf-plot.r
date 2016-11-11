# Performance Plots
# 
# Plot the performance of a GBMFit object produced by calling
# \code{\link{gbmt}}. This function is used in
# \code{\link{plot.GBMTPerformance}}.
# 
# 
# @param gbm_fit_obj a \code{GBMFit} created from an initial call to
# \code{\link{gbmt}}.
# 
# @param best_iter iteration specifying the optimum number of
# iterations. This is determined in \code{\link{gbmt_performance}}.
# 
# @param out_of_bag_curve logical indicating whether to plot the
# out-of-bag performance measures in a second plot.
# 
# @param overlay if TRUE and out_of_bag_curve=TRUE then a right y-axis
# is added to the training and test error plot and the estimated
# cumulative improvement in the loss function is plotted versus the
# iteration number.
# 
# @param method indicate the method used to estimate the optimal
# number of boosting iterations. \code{method="OOB"} computes the
# out-of-bag estimate and \code{method="test"} uses the test (or
# validation) dataset to compute an out-of-sample
# estimate. \code{method="cv"} extracts the optimal number of
# iterations using cross-validation if \code{gbmt} was called with
# \code{cv_folds}>1.
# 
# @param main the main title for the plot.
perf_plot <- function(gbm_fit_obj, best_iter, out_of_bag_curve,
                      overlay, method, main) {
  # Check inputs
  check_if_gbm_fit(gbm_fit_obj)
  if(!is.logical(overlay) || (length(overlay)) > 1 || is.na(overlay))
    stop("overlay must be a logical - excluding NA")
  
  if(!is.logical(out_of_bag_curve) || (length(out_of_bag_curve)) > 1 || is.na(out_of_bag_curve))
    stop("out_of_bag_curve must be a logical - excluding NA")
  
  par(mar=c(5,4,4,4)+.1)
  
  # Get y-axis label and limits
  ylab <- get_ylabel(gbm_fit_obj$distribution)
  if(gbm_fit_obj$params$train_fraction==1) {
    ylim <- switch(method,
                   cv=range(iteration_error(gbm_fit_obj, 'train'),
                       iteration_error(gbm_fit_obj, 'cv')),
                   test=range(iteration_error(gbm_fit_obj, 'train'),
                       iteration_error(gbm_fit_obj, 'valid')),
                   OOB=range(iteration_error('train')))
  } else {
        ylim <- range(iteration_error(gbm_fit_obj, 'train'),
                      iteration_error(gbm_fit_obj, 'valid'))
  }
  
  # Initial plot
  plot(iteration_error(gbm_fit_obj, 'train'),
       ylim=ylim,
       type="l",
       xlab="Iteration", ylab=ylab, main=main)
  
  if(gbm_fit_obj$params$train_fraction != 1) {
    lines(iteration_error(gbm_fit_obj, 'valid'), col="red")
  }
  if(method=="cv") {
    lines(iteration_error(gbm_fit_obj, 'cv'), col="green")
  }
  if(!is.na(best_iter)) abline(v=best_iter,col="blue",lwd=2,lty=2)
  
  # Plot out of bag curve
  if(out_of_bag_curve)
    plot_oobag(gbm_fit_obj, best_iter, overlay, ylab)
}
########## HELPER FUNCTIONS FOR PLOTTING ############
get_ylabel <- function(distribution) {
  UseMethod("get_ylabel", distribution)
}

get_ylabel.default <- function(distribution) {
  stop("distribution object not recognised - cannot get y label for plot")
}

get_ylabel.AdaBoostGBMDist <-function(distribution) {
  return("AdaBoost exponential bound")
}
get_ylabel.BernoulliGBMDist <- function(distribution) {
  return("Bernoulli deviance")
}
get_ylabel.CoxPHGBMDist <- function(distribution) {
  return("Cox partial deviance")
}
get_ylabel.GammaGBMDist <- function(distribution) {
  return("Gamma deviance")
}
get_ylabel.GaussianGBMDist <- function(distribution) {
  return("Squared error loss")
}
get_ylabel.HuberizedGBMDist <- function(distribution) {
  return("Hinged loss")
}
get_ylabel.LaplaceGBMDist <- function(distribution) {
  return("Absolute loss")
}
get_ylabel.PairwiseGBMDist <- function(distribution) {
  ylab <- switch(distribution$metric,
                 conc ="Fraction of concordant pairs",
                 ndcg="Normalized discounted cumulative gain",
                 map ="Mean average precision",
                 mrr ="Mean reciprocal rank")
  return(ylab)
}
get_ylabel.PoissonGBMDist <- function(distribution) {
  return("Poisson deviance")
}
get_ylabel.QuantileGBMDist <- function(distribution) {
  return("Quantile loss")
}
get_ylabel.TDistGBMDist <- function(distribution) {
  return("t-distribution deviance")
}
get_ylabel.TweedieGBMDist <- function(distribution) {
  return("Tweedie deviance")
}


plot_oobag <- function(gbm_fit_obj, best_iter, overlay, ylab) {
  # Get smoother
  smoother <- generate_smoother_oobag(gbm_fit_obj)
  
  # Plot smoothed out of bag improvement
  if(overlay) {
    par(new=TRUE)
    plot(smoother$x,
         cumsum(smoother$y),
         col="blue",
         type="l",
         xlab="",ylab="",
         axes=FALSE)
    axis(4,srt=0)
    at <- mean(range(smoother$y))
    mtext(paste("OOB improvement in",ylab),side=4,srt=270,line=2)
    abline(h=0,col="blue",lwd=2)
  }
  
  # Plot original out of bag improvement
  plot(gbm_fit_obj$oobag.improve,type="l",
       xlab="Iteration",
       ylab=paste("OOB change in",ylab))
  lines(smoother,col="red",lwd=2)
  abline(h=0,col="blue",lwd=1)
  abline(v=best_iter,col="blue",lwd=1)
}
