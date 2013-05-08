\name{gbm.object}
\alias{gbm.object}
\title{Generalized Boosted Regression Model Object}
\description{These are objects representing fitted \code{gbm}s.}
\section{Structure}{The following components must be included in a legitimate \code{gbm} object.}
\value{
\item{initF}{the "intercept" term, the initial predicted value to which trees
make adjustments}
\item{fit}{a vector containing the fitted values on the scale of regression
function (e.g. log-odds scale for bernoulli, log scale for poisson)}
\item{train.error}{a vector of length equal to the number of fitted trees
containing the value of the loss function for each boosting iteration
evaluated on the training data}
\item{valid.error}{a vector of length equal to the number of fitted trees
containing the value of the loss function for each boosting iteration
evaluated on the validation data}
\item{cv.error}{if \code{cv.folds}<2 this component is NULL. Otherwise, this 
component is a vector of length equal to the number of fitted trees
containing a cross-validated estimate of the loss function for each boosting 
iteration}
\item{oobag.improve}{a vector of length equal to the number of fitted trees
containing an out-of-bag estimate of the marginal reduction in the expected
value of the loss function. The out-of-bag estimate uses only the training
data and is useful for estimating the optimal number of boosting iterations.
See \code{\link{gbm.perf}}}
\item{trees}{a list containing the tree structures. The components are best
viewed using \code{\link{pretty.gbm.tree}}}
\item{c.splits}{a list of all the categorical splits in the collection of
trees. If the \code{trees[[i]]} component of a \code{gbm} object describes a
categorical split then the splitting value will refer to a component of
\code{c.splits}. That component of \code{c.splits} will be a vector of length
equal to the number of levels in the categorical split variable. -1 indicates
left, +1 indicates right, and 0 indicates that the level was not present in the
training data}
\item{cv.fitted}{If cross-validation was performed, the cross-validation
  predicted values on the scale of the linear predictor. That is, the
  fitted values from the ith CV-fold, for the model having been trained on
  the data in all other folds.}
}
\author{Greg Ridgeway \email{gregridgeway@gmail.com}}
\seealso{
\code{\link{gbm}}
}

\keyword{methods}
