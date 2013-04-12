\name{print.gbm}
\alias{print.gbm}
\alias{show.gbm}

\title{ Print model summary }
\description{
  Display basic information about a \code{gbm} object.
}
\usage{
\method{print}{gbm}(x, ...)
\method{show}{gbm}(x, ...)
}
\arguments{
  \item{x}{ an object of class \code{gbm}. }
  \item{\dots}{ arguments passed to \code{print.default}. }
}
\details{
  Prints some information about the model object. In particular, this method
  prints the call to \code{gbm()}, the type of loss function
  that was used, and the total number of iterations.

  If cross-validation was performed, the 'best' number of trees as
  estimated by cross-validation error is displayed. If a test set
  was used, the 'best' number
  of trees as estimated by the test set error is displayed.

  The number of available predictors, and the number of those having
  non-zero influence on predictions is given (which might be interesting
  in data mining applications).

  If multinomial, bernoulli or adaboost was used,
  the confusion matrix and prediction accuracy are printed (objects
  being allocated to the class with highest probability for multinomial
  and bernoulli). These classifications are performed on the entire
  training
  data using the model with the 'best' number of trees as described
  above, or the maximum number of trees if the 'best' cannot be
  computed.

  If the 'distribution' was specified as gaussian, laplace, quantile
  or t-distribution, a summary of the residuals is displayed.
  The residuals are for the training data with the model at the 'best'
  number of trees, as
  described above, or the maximum number of trees if the 'best' cannot
  be computed.
}


\author{ Harry Southworth, Daniel Edwards }

\seealso{ \code{\link{gbm}} }
\examples{
data(iris)
iris.mod <- gbm(Species ~ ., distribution="multinomial", data=iris,
                 n.trees=2000, shrinkage=0.01, cv.folds=5,
                 verbose=FALSE, n.cores=1)
iris.mod
#data(lung)
#lung.mod <- gbm(Surv(time, status) ~ ., distribution="coxph", data=lung,
#                 n.trees=2000, shrinkage=0.01, cv.folds=5,verbose =FALSE)
#lung.mod
}
\keyword{models}
\keyword{nonlinear}
\keyword{survival}
\keyword{nonparametric}
