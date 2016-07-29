#' Baseline hazard function
#' 
#' Computes the Breslow estimator of the baseline hazard function for a
#' proportional hazard regression model - only for censored survival data.
#' 
#' The proportional hazard model assumes h(t|x)=lambda(t)*exp(f(x)).
#' \code{\link{gbm2}} can estimate the f(x) component via partial likelihood.
#' After estimating f(x), \code{baseline_hazard} can compute a nonparametric
#' estimate of lambda(t).
#' 
#' @usage baseline_hazard(surv_times, delta, coxph_pred, eval_times=NULL, smooth=FALSE,
#' cumulative=TRUE)
#' 
#' @param surv_times the survival times - an atomic vector of doubles
#' 
#' @param delta the censoring indicator - a vector same length as surv_times
#' 
#' @param coxph_pred the predicted values of the regression model on the log hazard
#' scale
#' 
#' @param eval_times values at which the baseline hazard will be evaluated
#' 
#' @param smooth if \code{TRUE} \code{baseline_hazard} will smooth the estimated
#' baseline hazard using Friedman's super smoother \code{\link{supsmu}}
#' 
#' @param cumulative if \code{TRUE} the cumulative survival function will be
#' computed
#' 
#' @return a vector of length equal to the length of surv_times (or of length
#' \code{eval_times} if \code{eval_times} is not \code{NULL}) containing the baseline
#' hazard evaluated at t (or at \code{eval_times} if \code{eval_times} is not
#' \code{NULL}). If \code{cumulative} is set to \code{TRUE} then the returned
#' vector evaluates the cumulative hazard function at those values.
#' 
#' @author James Hickey, Greg Ridgeway \email{gregridgeway@@gmail.com}
#' 
#' @seealso \code{\link[survival]{survfit}}, \code{\link{gbm2}}
#' 
#' @references N. Breslow (1972). "Discussion of `Regression Models and
#' Life-Tables' by D.R. Cox," Journal of the Royal Statistical Society, Series
#' B, 34(2):216-217.
#' 
#' N. Breslow (1974). "Covariance analysis of censored survival data,"
#' Biometrics 30:89-99.
#' 
#' @keywords methods survival
#' 
#' @export

baseline_hazard <- function(surv_times, delta, coxph_preds,
                            eval_times=NULL,
                            smooth=FALSE,
                            cumulative=TRUE) {
  # Initial checks
  if(!is.logical(cumulative) || (length(cumulative) > 1)) {
    stop("cumulative must be a logical")
  }
  if(!is.logical(smooth) || (length(smooth) > 1)) {
    stop("smooth must be a logical")
  }
  if(!is.atomic(surv_times) || !any(surv_times == as.double(surv_times))) {
    stop("surv_times must be a vector of doubles")
  }
  
  if(!is.atomic(delta) || any(!(delta %in% c(0, 1))) ) {
    stop("delta must be a vector with elements eiter 0 or 1")
  }
  
  if(!is.null(eval_times) && 
     (!is.atomic(eval_times) || !any(eval_times == as.double(eval_times))) ) {
    stop("eval_times must be a vector of doubles if not NULL")
  }
  
  if(!smooth && !cumulative)
  {
    if(!is.null(eval_times))
    {
      stop("Cannot evaluate unsmoothed baseline hazard at eval_times.")
    }
  }
  
  # Calculate survival function - alpha
  unique_death_times <- sort(unique(surv_times[delta==1]))
  alpha <- length(unique_death_times)
  for(i in seq_len(length(unique_death_times))) {
    alpha[i] <- sum(surv_times[delta==1]==unique_death_times[i])/
      sum(exp(coxph_preds[surv_times >= unique_death_times[i]]))
  }
  
  # Check if need cumulative survival function
  if(cumulative) alpha <- cumsum(alpha)
  
  # Set evaluation times
  if(!is.null(eval_times)) eval_times <- surv_times
  
  # Use super smoother and evaluate 
  lambda <- list(x=unique_death_times, y=alpha)
  if(smooth) {
    lambda <- supsmu(unique_death_times, alpha)
  } 

  return(approx(lambda$x, lambda$y, xout=eval_times)$y)
}