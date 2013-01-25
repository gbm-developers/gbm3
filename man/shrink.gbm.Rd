\name{shrink.gbm}
\alias{shrink.gbm}
\title{ L1 shrinkage of the predictor variables in a GBM }
\description{
Performs recursive shrinkage in each of the trees in a GBM fit using different shrinkage parameters for each variable.
}
\usage{
shrink.gbm(object, 
           n.trees, 
           lambda = rep(10, length(object$var.names)), 
           ...)
}
\arguments{
  \item{object}{ A \code{\link{gbm.object}} }
  \item{n.trees}{ the number of trees to use }
  \item{lambda}{ a vector with length equal to the number of variables containing the shrinkage parameter for each variable }
  \item{\dots}{ other parameters (ignored) }
}
\details{
This function is currently experimental. Used in conjunction with a gradient ascent search for inclusion of variables.
}
\value{
  \item{predF}{Predicted values from the shrunken tree}
  \item{objective}{The value of the loss function associated with the predicted values}
  \item{gradient}{A vector with length equal to the number of variables containing the derivative of the objective function with respect to beta, the logit transform of the shrinkage parameter for each variable}
}
\references{ Hastie, T. J., and Pregibon, D. \href{http://www-stat.stanford.edu/~hastie/Papers/shrinktree.ps}{Shrinking Trees}. AT&T Bell Laboratories Technical Report (March 1990).}
\author{ Greg Ridgeway \email{gregridgeway@gmail.com} }


\section{Warning}{This function is experimental.}

\seealso{ \code{\link{shrink.gbm.pred}}, \code{\link{gbm}} }
\keyword{ methods}% at least one, from doc/KEYWORDS

