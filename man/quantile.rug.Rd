\name{quantile.rug}
\alias{quantile.rug}
\title{Quantile rug plot}
\description{Marks the quantiles on the axes of the current plot.}
\usage{
quantile.rug(x,prob=(0:10)/10,...)
}
\arguments{
\item{x}{a numeric vector.}
\item{prob}{the quantiles of x to mark on the x-axis.}
\item{...}{additional graphics parameters currently ignored.}
}
\value{No return values}
\author{Greg Ridgeway \email{gregridgeway@gmail.com}}
\seealso{
\code{\link[graphics]{plot}}, 
\code{\link[stats]{quantile}}, 
\code{\link[base]{jitter}}, 
\code{\link[graphics]{rug}}.
}
\examples{
x <- rnorm(100)
y <- rnorm(100)
plot(x,y)
quantile.rug(x)
}
\keyword{aplot}
