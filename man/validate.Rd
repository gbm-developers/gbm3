\name{validate.gbm}
\alias{validate.gbm}
\alias{test.gbm}
\alias{test.relative.influence}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Test the \code{gbm} package.}
\description{
Run tests on \code{gbm} functions to perform logical checks and
reproducibility.
}
\usage{
validate.gbm()
}
%- maybe also 'usage' for other objects documented here.
\details{
The function uses functionality in the \code{RUnit} package.
A fairly small validation suite is executed that checks to see that
relative influence identifies sensible variables from simulated data,
and that predictions from GBMs with Gaussian, Cox or binomial distributions
are sensible,
}
\value{
An object of class \code{RUnitTestData}. See the help for \code{RUnit} for
details.
}
\author{
Harry Southworth
}
\note{The test suite is not comprehensive.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{gbm}}
}
\examples{
# Uncomment the following lines to run - commented out to make CRAN happy
#library(RUnit)
#val <- validate.texmex()
#printHTMLProtocol(val, "texmexReport.html")
}
\keyword{models}

