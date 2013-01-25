\name{gbm.roc.area}
\alias{gbm.roc.area}
\alias{gbm.conc}
\alias{ir.measure.conc}
\alias{ir.measure.auc}
\alias{ir.measure.mrr}
\alias{ir.measure.map}
\alias{ir.measure.ndcg}
\alias{perf.pairwise}
\title{
Compute Information Retrieval measures.
}
\description{
 Functions to compute Information Retrieval measures for pairwise loss for
 a single group. The function returns the respective metric, or a negative value if
it is undefined for the given group. 
}
\usage{
gbm.roc.area(obs, pred)
ir.measure.conc(y.f, max.rank)
ir.measure.auc(y.f, max.rank)
ir.measure.mrr(y.f, max.rank)
ir.measure.map(y.f, max.rank)
ir.measure.ndcg(y.f, max.rank)
perf.pairwise(y, f, group, metric="ndcg", w=NULL, max.rank=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs}{Observed value}
  \item{pred}{Predicted value}
\item{metric}{What type of performance measure to compute.}
\item{y, y.f, f, w, group, max.rank}{Used internally.}
}
\details{
For simplicity, we have no special handling for ties;
instead, we break ties randomly. This is slightly
inaccurate for individual groups, but should have
only a small effect on the overall measure.

\code{gbm.conc} computes the concordance index:
Fraction of all pairs (i,j) with i<j, x[i] != x[j], such that x[j] < x[i]

If \code{obs} is binary, then
\code{gbm.roc.area(obs, pred) = gbm.conc(obs[order(-pred)])}.

\code{gbm.conc} is more general as it allows non-binary targets,
but is significantly slower.
}
\value{The requested performance measure.}
\references{
C. Burges (2010). "From RankNet to LambdaRank to LambdaMART: An Overview",
Microsoft Research Technical Report MSR-TR-2010-82.
}
\author{
Stefan Schroedl
}

\seealso{
\code{\link{gbm}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
\keyword{ models }

