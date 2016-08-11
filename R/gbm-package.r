#' @importFrom Rcpp sourceCpp
#' @useDynLib gbm
#' @importFrom survival Surv
#' @import lattice
#' @import parallel
#' @importFrom graphics abline axis barplot lines mtext par plot polygon rug
#' @importFrom grDevices rainbow
#' @importFrom stats approx binomial delete.response gaussian glm loess
#' @importFrom stats mad median model.extract model.frame model.offset
#' @importFrom stats model.response model.weights na.pass poisson
#' @importFrom stats predict quantile reformulate runif supsmu
#' @importFrom stats terms weighted.mean
NULL

## on package load set default options
.onLoad <- function(libname, pkgname) {
    already.set <- options()

    if (!('gbm.parallel' %in% names(already.set))) {
        options(gbm.parallel=gbmParallel())
    }

    invisible()
}

#' gbm internal functions
#' 
#' Helper functions for preprocessing data prior to building the model
#' 
#' These are functions used internally by \code{gbmt} and not intended for
#' direct use by the user.
#' 
#' @name gbm-internal
NULL





#' Generalized Boosted Regression Model Object
#' 
#' These are objects representing fitted \code{gbm}s.
#' 
#' 
#' @return \item{initF}{the "intercept" term, the initial predicted value to
#' which trees make adjustments} 
#' \item{fit}{a vector containing the fitted
#' values on the scale of regression function (e.g. log-odds scale for
#' bernoulli, log scale for poisson)} 
#' \item{train.error}{a vector of length
#' equal to the number of fitted trees containing the value of the loss
#' function for each boosting iteration evaluated on the training data}
#' \item{valid.error}{a vector of length equal to the number of fitted trees
#' containing the value of the loss function for each boosting iteration
#' evaluated on the validation data}
#' \item{cv_error}{if \code{cv_folds}<2 this
#' component is NULL. Otherwise, this component is a vector of length equal to
#' the number of fitted trees containing a cross-validated estimate of the loss
#' function for each boosting iteration}
#' \item{oobag.improve}{a vector of
#' length equal to the number of fitted trees containing an out-of-bag estimate
#' of the marginal reduction in the expected value of the loss function. The
#' out-of-bag estimate uses only the training data and is useful for estimating
#' the optimal number of boosting iterations. See \code{\link{gbm_perf}}}
#' \item{trees}{a list containing the tree structures. The components are best
#' viewed using \code{\link{pretty_gbm_tree}}}
#' \item{c.splits}{a list of all
#' the categorical splits in the collection of trees. If the \code{trees[[i]]}
#' component of a \code{gbm} object describes a categorical split then the
#' splitting value will refer to a component of \code{c.splits}. That component
#' of \code{c.splits} will be a vector of length equal to the number of levels
#' in the categorical split variable. -1 indicates left, +1 indicates right,
#' and 0 indicates that the level was not present in the training data}
#' \item{cv_fitted}{If cross-validation was performed, the cross-validation
#' predicted values on the scale of the linear predictor. That is, the fitted
#' values from the ith CV-fold, for the model having been trained on the data
#' in all other folds.}
#' @section Structure: The following components must be included in a
#' legitimate \code{GBMFit} object.
#' 
#' @author Greg Ridgeway \email{gregridgeway@@gmail.com}
#' @seealso \code{\link{gbmt}}
#' @keywords methods
#' @name gbm_object
NULL


#' Generalized Boosted Regression Models
#' 
#' This package implements extensions to Freund and Schapire's AdaBoost
#' algorithm and J. Friedman's gradient boosting machine. Includes regression
#' methods for least squares, absolute loss, logistic, Poisson, Cox
#' proportional hazards partial likelihood, t-distribution,
#' AdaBoost exponential loss, Learning to Rank, and Huberized hinge loss.
#' 
#' \tabular{ll}{Package: \tab gbm\cr Version: \tab 2.1-0.6\cr Date: \tab
#' 2014-08-12\cr Depends: \tab R (>= 2.9.0), survival, lattice, mgcv\cr
#' License: \tab GPL (version 2 or newer)\cr} Index:
#' \preformatted{baseline_hazard Baseline hazard function calibrate_plot
#' Calibration plot gbmt Generalized Boosted Regression Modeling gbm_perf GBM performance
#' gbmt_plot Marginal plots of fitted gbm objects predict.GBMFit Predict method for
#' GBM Model Fits pretty_gbm_tree Print gbm tree components quantile_rug
#' Quantile rug plot relative_influence Methods for estimating relative
#' influence}
#' 
#' Further information is available in the following vignettes: \tabular{ll}{
#' \code{gbm} \tab Generalized Boosted Models: A guide to the gbm package
#' (source, pdf)\cr}
#' 
#' @name gbm-package
#' @docType package
#' @author James Hickey, Greg Ridgeway \email{gregridgeway@@gmail.com} with contributions by
#' Daniel Edwards, Brian Kriegler, Stefan Schroedl and Harry Southworth.
#' @references Y. Freund and R.E. Schapire (1997) \dQuote{A decision-theoretic
#' generalization of on-line learning and an application to boosting,}
#' \emph{Journal of Computer and System Sciences,} 55(1):119-139.
#' 
#' G. Ridgeway (1999). \dQuote{The state of boosting,} \emph{Computing Science
#' and Statistics} 31:172-181.
#' 
#' J.H. Friedman, T. Hastie, R. Tibshirani (2000). \dQuote{Additive Logistic
#' Regression: a Statistical View of Boosting,} \emph{Annals of Statistics}
#' 28(2):337-374.
#' 
#' J.H. Friedman (2001). \dQuote{Greedy Function Approximation: A Gradient
#' Boosting Machine,} \emph{Annals of Statistics} 29(5):1189-1232.
#' 
#' J.H. Friedman (2002). \dQuote{Stochastic Gradient Boosting,}
#' \emph{Computational Statistics and Data Analysis} 38(4):367-378.
#' 
#' The \href{http://www-stat.stanford.edu/~jhf/R-MART.html}{MART} website.
#' @keywords package
NULL



