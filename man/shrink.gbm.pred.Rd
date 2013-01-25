\name{shrink.gbm.pred}
\alias{shrink.gbm.pred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Predictions from a shrunked GBM }
\description{
  Makes predictions from a shrunken GBM model.
}
\usage{
shrink.gbm.pred(object, 
                newdata, 
                n.trees, 
                lambda = rep(1, length(object$var.names)), 
                ...)
}

\arguments{
  \item{object}{ a \code{\link{gbm.object}} }
  \item{newdata}{ dataset for predictions }
  \item{n.trees}{ the number of trees to use }
  \item{lambda}{ a vector with length equal to the number of variables containing the shrinkage parameter for each variable }
  \item{\dots}{ other parameters (ignored) }
}

\value{
A vector with length equal to the number of observations in newdata containing the predictions
}

\author{ Greg Ridgeway \email{gregridgeway@gmail.com} }

\section{Warning}{This function is experimental}

\seealso{ \code{\link{shrink.gbm}}, \code{\link{gbm}} }

\keyword{ methods }
