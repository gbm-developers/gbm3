\name{basehaz.gbm}
\alias{basehaz.gbm}
\title{ Baseline hazard function }
\description{
  Computes the Breslow estimator of the baseline hazard function for a proportional hazard regression model
}
\usage{
basehaz.gbm(t, delta, f.x, 
            t.eval = NULL, 
            smooth = FALSE, 
            cumulative = TRUE)
}
\arguments{
  \item{t}{ the survival times }
  \item{delta}{ the censoring indicator }
  \item{f.x}{ the predicted values of the regression model on the log hazard scale }
  \item{t.eval}{ values at which the baseline hazard will be evaluated }
  \item{smooth}{ if \code{TRUE} \code{basehaz.gbm} will smooth the estimated baseline hazard using Friedman's super smoother \code{\link{supsmu}}}
  \item{cumulative}{ if \code{TRUE} the cumulative survival function will be computed }
}
\details{
  The proportional hazard model assumes h(t|x)=lambda(t)*exp(f(x)). \code{\link{gbm}} can estimate the f(x) component via partial likelihood. After estimating f(x), \code{basehaz.gbm} can compute the a nonparametric estimate of lambda(t).
}
\value{
 a vector of length equal to the length of t (or of length \code{t.eval} if \code{t.eval} is not \code{NULL}) containing the baseline hazard evaluated at t (or at \code{t.eval} if \code{t.eval} is not \code{NULL}). If \code{cumulative} is set to \code{TRUE} then the returned vector evaluates the cumulative hazard function at those values.
}
\references{N. Breslow (1972). "Disussion of `Regression Models and Life-Tables' by D.R. Cox," Journal of the Royal Statistical Society, Series B, 34(2):216-217.

N. Breslow (1974). "Covariance analysis of censored survival data," Biometrics 30:89-99.
} 
\author{ Greg Ridgeway \email{gregridgeway@gmail.com}}



\seealso{ \code{\link[survival]{survfit}}, \code{\link{gbm}} }

\keyword{ methods }
\keyword{ survival }
