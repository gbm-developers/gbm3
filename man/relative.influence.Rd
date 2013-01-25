\name{relative.influence}
\alias{relative.influence}
\alias{permutation.test.gbm}
\alias{gbm.loss}
\title{ Methods for estimating relative influence }
\description{
Helper functions for computing the relative influence of each variable in the gbm object.
}
\usage{
relative.influence(object, n.trees, scale., sort.)
permutation.test.gbm(object, n.trees)
gbm.loss(y,f,w,offset,dist,baseline, group, max.rank)
}
\arguments{
\item{object}{a \code{gbm} object created from an initial call to \code{\link{gbm}}.}
\item{n.trees}{ the number of trees to use for computations. If not provided, the
    the function will guess: if a test set was used in fitting, the number of
    trees resulting in lowest test set error will be used; otherwise, if
    cross-validation was performed, the number of trees resulting in lowest
    cross-validation error will be used; otherwise, all trees will be used.}
\item{scale.}{ whether or not the result should be scaled. Defaults to \code{FALSE}.}
\item{sort.}{ whether or not the results should be (reverse) sorted.
         Defaults to \code{FALSE}.}
\item{y,f,w,offset,dist,baseline}{For \code{gbm.loss}: These components are the
outcome, predicted value, observation weight, offset, distribution, and comparison
loss function, respectively.}
\item{group, max.rank}{Used internally when \code{distribution = \'pairwise\'}.}
}
\details{
This is not intended for end-user use. These functions offer the different
methods for computing the relative influence in \code{\link{summary.gbm}}.
\code{gbm.loss} is a helper function for \code{permutation.test.gbm}.
}
\value{
By default, returns an unprocessed vector of estimated relative influences.
If the \code{scale.} and \code{sort.} arguments are used, returns a processed
version of the same.
}
\references{
J.H. Friedman (2001). "Greedy Function Approximation: A Gradient Boosting
Machine," Annals of Statistics 29(5):1189-1232.

L. Breiman (2001). \href{http://oz.berkeley.edu/users/breiman/randomforest2001.pdf}{Random Forests}.

}
\author{Greg Ridgeway \email{gregridgeway@gmail.com}}

\seealso{ \code{\link{summary.gbm}} }

\keyword{ hplot }
