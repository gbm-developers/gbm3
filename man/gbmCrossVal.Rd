\name{gbmCrossVal}
\alias{gbmCrossVal}
\alias{gbmCrossValModelBuild}
\alias{gbmDoFold}
\alias{gbmCrossValErr}
\alias{gbmCrossValPredictions}
\title{Cross-validate a gbm}
\description{Functions for cross-validating gbm. These functions are
  used internally and are not intended for end-user direct usage.}
\usage{
gbmCrossVal(cv.folds, nTrain, n.cores,
class.stratify.cv, data,
x, y, offset, distribution, w, var.monotone,
n.trees, interaction.depth, n.minobsinnode,
shrinkage, bag.fraction,
var.names, response.name, group)

gbmCrossValModelBuild(cv.folds, cv.group, n.cores,
i.train, x, y, offset,
distribution, w, var.monotone,
n.trees, interaction.depth,
n.minobsinnode, shrinkage,
bag.fraction, var.names,
response.name, group)

gbmDoFold(X, i.train, x, y, offset, distribution, w, var.monotone, n.trees,
         interaction.depth, n.minobsinnode, shrinkage, bag.fraction,
         cv.group, var.names, response.name, group, s)

gbmCrossValErr(cv.models, cv.folds, cv.group, nTrain, n.trees)

gbmCrossValPredictions(cv.models, cv.folds, cv.group,
best.iter.cv, distribution, data, y)


}
\arguments{
\item{cv.folds}{The number of cross-validation folds.}
\item{nTrain}{The number of training samples.}
\item{n.cores}{The number of cores to use.}
\item{class.stratify.cv}{Whether or not stratified cross-validation
      samples are used.}
\item{data}{The data.}
\item{x}{The model matrix.}
\item{y}{The response variable.}
\item{offset}{The offset.}
\item{distribution}{The type of loss function. See \code{\link{gbm}}.}
\item{w}{Observation weights.}
\item{var.monotone}{See \code{\link{gbm}}.}
\item{n.trees}{The number of trees to fit.}
\item{interaction.depth}{The degree of allowed interactions. See \code{\link{gbm}}.}
\item{n.minobsinnode}{See \code{\link{gbm}}.}
\item{shrinkage}{See \code{\link{gbm}}.}
\item{bag.fraction}{See \code{\link{gbm}}.}
\item{var.names}{See \code{\link{gbm}}.}
\item{response.name}{See \code{\link{gbm}}.}
\item{group}{Used when \code{distribution = "pairwise"}. See
  \code{\link{gbm}}.}
\item{i.train}{Items in the training set.}
\item{cv.models}{A list containing the models for each fold.}
\item{cv.group}{A vector indicating the cross-validation fold for each
  member of the training set.}
\item{best.iter.cv}{The iteration with lowest cross-validation error.}
\item{X}{Index (cross-validation fold) on which to subset.}
\item{s}{Random seed.}

} % Close arguments
\details{
These functions are not intended for end-user direct usage, but are used
internally by \code{gbm}.}
\value{A list containing the cross-validation error and predictions.}
\references{
J.H. Friedman (2001). "Greedy Function Approximation: A Gradient Boosting
Machine," Annals of Statistics 29(5):1189-1232.

L. Breiman (2001). \href{http://oz.berkeley.edu/users/breiman/randomforest2001.pdf}{Random Forests}.

}
\author{Greg Ridgeway \email{gregridgeway@gmail.com}}

\seealso{ \code{\link{gbm}} }

\keyword{ models }
