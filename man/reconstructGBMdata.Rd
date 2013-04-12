\name{reconstructGBMdata}
\Rdversion{1.1}
\alias{reconstructGBMdata}
\title{Reconstruct a GBM's Source Data}


\description{
Helper function to reconstitute the data for plots and summaries. This
function is not intended for the user to call directly.
}
\usage{
reconstructGBMdata(x)
}
\arguments{
  \item{x}{
a \code{\link{gbm.object}} initially fit using \code{\link{gbm}}
}
}
\value{
Returns a data used to fit the gbm in a format that can subsequently be used
for plots and summaries
}
\author{
Harry Southworth
}

\seealso{
\code{\link{gbm}}, \code{\link{gbm.object}}
}

\keyword{ manip }
