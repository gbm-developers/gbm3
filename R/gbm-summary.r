#' Summary of a GBMFit object
#' 
#' Computes the relative influence of each variable in the \code{GBMFit} object.
#' 
#' For \code{GBMGaussianDist} this returns exactly the reduction of
#' squared error attributable to each variable. For other loss functions this
#' returns the reduction attributable to each variable in sum of squared error
#' in predicting the gradient on each iteration. It describes the relative
#' influence of each variable in reducing the loss function. See the references
#' below for exact details on the computation.
#' 
#' @usage summary(gbm_fit_obj, cBars=length(gbm_fit_obj$variables$var_names), 
#'                num_trees=length(gbm_fit_obj$trees), plot_it=TRUE, order_it=TRUE,
#'                method=relative_influence, normalize=TRUE, ...)
#' 
#' @param gbm_fit_obj a \code{GBMFit} object created from an initial call to
#' \code{\link{gbm2}}.
#' 
#' @param cBars the number of bars to plot. If \code{order_it=TRUE} then only the
#' \code{cBars} variables with the largest relative influence will appear in
#' the barplot. If \code{order_it=FALSE} then the first \code{cBars} variables
#' will appear in the plot. In either case, the function will return the
#' relative influence of all of the variables.
#' 
#' @param num_trees the number of trees used to generate the plot. Only the first
#' \code{num_trees} trees will be used.
#' 
#' @param plot_it an indicator as to whether the plot is generated.
#' 
#' @param order_it an indicator as to whether the plotted and/or returned relative
#' influences are sorted.
#' 
#' @param method The function used to compute the relative influence.
#' \code{\link{relative_influence}} is the default and is the same as that
#' described in Friedman (2001). The other current (and experimental) choice is
#' \code{\link{permutation_test_gbm}}. This method randomly permutes each
#' predictor variable at a time and computes the associated reduction in
#' predictive performance. This is similar to the variable importance measures
#' Breiman uses for random forests, but \code{gbm} currently computes using the
#' entire training dataset (not the out-of-bag observations).
#' 
#' @param normalize if \code{FALSE} then \code{summary.gbm} returns the
#' unnormalized influence.
#' 
#' @param ...  other arguments passed to the plot function.
#' 
#' @return Returns a data frame where the first component is the variable name
#' and the second is the computed relative influence, normalized to sum to 100.
#' @author Greg Ridgeway \email{gregridgeway@@gmail.com}
#' @seealso \code{\link{gbm2}}
#' @references J.H. Friedman (2001). "Greedy Function Approximation: A Gradient
#' Boosting Machine," Annals of Statistics 29(5):1189-1232.
#' 
#' L. Breiman (2001). \href{http://oz.berkeley.edu/users/breiman/randomforest2001.pdf}{Random Forests}.
#' @keywords hplot
#' @export
#' 

summary.GBMFit <- function(gbm_fit_obj,
                        cBars=length(gbm_fit_obj$variables$var_names),
                        num_trees=length(gbm_fit_obj$trees),
                        plot_it=TRUE,
                        order_it=TRUE,
                        method=relative_influence,
                        normalize=TRUE,
                        ...)
{
  # Initial checks
  check_if_natural_number(num_trees)
  check_if_natural_number(cBars)
  check_if_gbm_fit(gbm_fit_obj)
  if(!is.logical(plot_it) || (length(plot_it) > 1)) {
    stop("argument plot_it must be a logical")
  }  
  if(!is.logical(order_it) || (length(order_it) > 1)) {
    stop("argument order_it must be a logical")
  }  
  if(!is.logical(normalize) || (length(normalize) > 1)) {
    stop("argument normalize must be a logical")
  }  
  
  # Set inputs (if required)
  if(cBars==0) cBars <- min(10, length(gbm_fit_obj$variables$var_names))
  if(cBars>length(gbm_fit_obj$variables$var_names)) cBars <- length(gbm_fit_obj$variables$var_names)
  if(num_trees > gbm_fit_obj$params$num_trees)
    warning("Exceeded total number of GBM terms. Results use num_trees=", gbm_fit_obj$params$num_trees," terms.\n")
  num_trees <- min(num_trees, gbm_fit_obj$params$num_trees)
  
  # Calculate relative influence and order/normalize
  rel_inf <- method(gbm_fit_obj, num_trees=num_trees)
  rel_inf[rel_inf<0] <- 0
  if(normalize) rel_inf <- 100*rel_inf/sum(rel_inf)
  
  ordering <- seq_len(length(rel_inf))
  if(order_it) {
    ordering <- order(-rel_inf)
  }
  
  # Bar plot of relative influence
  if(plot_it) {
    barplot(rel_inf[ordering[cBars:1]],
            horiz=TRUE,
            col=rainbow(cBars,start=3/6,end=4/6),
            names=gbm_fit_obj$variables$var_names[ordering[cBars:1]],
            xlab="Relative influence",
            las=1,...)
  }
  return(data.frame(var=gbm_fit_obj$variables$var_names[ordering],
                    rel_inf=rel_inf[ordering]))
}