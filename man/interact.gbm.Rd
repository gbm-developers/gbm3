\name{interact.gbm}
\alias{interact.gbm}
\title{ Estimate the strength of interaction effects }
\description{ Computes Friedman's H-statistic to assess the strength of variable interactions. }
\usage{
interact.gbm(x,
             data,
             i.var = 1,
             n.trees = x$n.trees)
}
\arguments{
  \item{x}{ a \code{\link{gbm.object}} fitted using a call to \code{\link{gbm}}}
  \item{data}{ the dataset used to construct \code{x}. If the original dataset is
     large, a random subsample may be used to accelerate the computation in
     \code{interact.gbm}}
  \item{i.var}{a vector of indices or the names of the variables for compute
     the interaction effect. If using indices, the variables are indexed in the
     same order that they appear in the initial \code{gbm} formula.}
  \item{n.trees}{ the number of trees used to generate the plot. Only the first
     \code{n.trees} trees will be used}
}
\details{
\code{interact.gbm} computes Friedman's H-statistic to assess the relative
strength of interaction effects in non-linear models. H is on the scale of
[0-1] with higher values indicating larger interaction effects. To connect to
a more familiar measure, if \eqn{x_1} and \eqn{x_2} are uncorrelated covariates
with mean 0 and variance 1 and the model is of the form
\deqn{y=\beta_0+\beta_1x_1+\beta_2x_2+\beta_3x_3}
then
\deqn{H=\frac{\beta_3}{\sqrt{\beta_1^2+\beta_2^2+\beta_3^2}}}

Note that if the main effects are weak, the estimated H will be unstable. For
example, if (in the case of a two-way interaction) neither main effect is in
the selected model (relative influence is zero), the result will be 0/0. Also,
with weak main effects, rounding errors can result in values of H > 1 which are
not possible.
}
\value{
Returns the value of \eqn{H}.
}
\references{
J.H. Friedman and B.E. Popescu (2005). \dQuote{Predictive Learning via Rule
Ensembles.} Section 8.1
}
\author{Greg Ridgeway \email{gregridgeway@gmail.com}}

\seealso{ \code{\link{gbm}}, \code{\link{gbm.object}} }

\keyword{ methods }
