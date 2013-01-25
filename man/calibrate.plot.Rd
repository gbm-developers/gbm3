\name{calibrate.plot}
\alias{calibrate.plot}
\title{Calibration plot}
\description{
An experimental diagnostic tool that plots the fitted values versus the actual average values.
Currently developed for only \code{distribution="bernoulli"}.
}
\usage{
calibrate.plot(y,p,
               distribution="bernoulli",
               replace=TRUE,
               line.par=list(col="black"),
               shade.col="lightyellow",
               shade.density=NULL,
               rug.par=list(side=1),
               xlab="Predicted value",
               ylab="Observed average",
               xlim=NULL,ylim=NULL,
               knots=NULL,df=6,
               ...)
}
\arguments{
  \item{y}{ the outcome 0-1 variable }
  \item{p}{ the predictions estimating E(y|x) }
  \item{distribution}{the loss function used in creating \code{p}.
                      \code{bernoulli} and \code{poisson} are currently the
                      only special options. All others default to squared error
                      assuming \code{gaussian}}
  \item{replace}{ determines whether this plot will replace or overlay the current plot.
                  \code{replace=FALSE} is useful for comparing the calibration of several
                  methods}
  \item{line.par}{ graphics parameters for the line }
  \item{shade.col}{ color for shading the 2 SE region. \code{shade.col=NA} implies no 2 SE
                    region}
  \item{shade.density}{ the \code{density} parameter for \code{\link{polygon}}}
  \item{rug.par}{graphics parameters passed to \code{\link{rug}}}
  \item{xlab}{x-axis label corresponding to the predicted values}
  \item{ylab}{y-axis label corresponding to the observed average}
  \item{xlim,ylim}{x and y-axis limits. If not specified te function will select
                   limits}
  \item{knots,df}{these parameters are passed directly to 
                  \code{\link[splines]{ns}} for constructing a natural spline 
                  smoother for the calibration curve}
  \item{...}{ other graphics parameters passed on to the plot function }
}
\details{
Uses natural splines to estimate E(y|p). Well-calibrated predictions
imply that E(y|p) = p. The plot also includes a pointwise 95% confidence
band.
}
\value{
\code{calibrate.plot} returns no values.
}
\references{
J.F. Yates (1982). "External correspondence: decomposition of the mean
probability score," Organisational Behaviour and Human Performance 30:132-156.

D.J. Spiegelhalter (1986). "Probabilistic Prediction in Patient Management
and Clinical Trials," Statistics in Medicine 5:421-433.
}
\author{Greg Ridgeway \email{gregridgeway@gmail.com}}
\examples{
# Don't want R CMD check to think there is a dependency on rpart
# so comment out the example
#library(rpart)
#data(kyphosis)
#y <- as.numeric(kyphosis$Kyphosis)-1
#x <- kyphosis$Age
#glm1 <- glm(y~poly(x,2),family=binomial)
#p <- predict(glm1,type="response")
#calibrate.plot(y, p, xlim=c(0,0.6), ylim=c(0,0.6))
}
\keyword{ hplot }
