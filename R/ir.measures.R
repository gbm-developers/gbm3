# Functions to compute IR measures for pairwise loss for
# a single group
# Notes:
# * Inputs are passed as a 2-elemen (y,f) list, to
#   facilitate the 'by' iteration
# * Return the respective metric, or a negative value if
#   it is undefined for the given group
# * For simplicity, we have no special handling for ties;
#   instead, we break ties randomly. This is slightly
#   inaccurate for individual groups, but should have
#   a small effect on the overall measure.


# Area under ROC curve = ratio of correctly ranking pairs
gbm.roc.area <- function(obs, pred)
{
   n1 <- sum(obs)
   n <- length(obs)
   if (n==n1) { return(1) }
   # Fraction of concordant pairs
   # = sum_{pos}(rank-1) / #pairs with different labels
   # #pairs = n1 * (n-n1)
   return ((mean(rank(pred)[obs > 0]) - (n1 + 1)/2)/(n - n1))
}

# Concordance Index:
# Fraction of all pairs (i,j) with i<j, x[i] != x[j], such that x[j] < x[i]
# Invariant: if obs is binary, then
#      gbm.roc.area(obs, pred) = gbm.conc(obs[order(-pred)])
# gbm.conc is more general as it allows non-binary targets,
# but is significantly slower
gbm.conc <- function(x)
{
   lx <- length(x)
   return (sum(mapply(function(r) { sum(x[(r+1):lx]<x[r]) }, 1:(lx-1))))
}

ir.measure.conc <- function(y.f, max.rank=0)
{
   # Note: max.rank is meaningless for CONC

   y           <- y.f[[1]]
   f           <- y.f[[2]]

   tab         <- table(y)
   csum        <- cumsum(tab)
   total.pairs <- sum(tab * (csum - tab))

   if (total.pairs == 0)
   {
      return (-1.0)
   }
   else
   {
      return (gbm.conc(y[order(-f)]) / total.pairs)
   }
}

ir.measure.auc <- function(y.f, max.rank=0)
{
   # Note: max.rank is meaningless for AUC
   y       <- y.f[[1]]
   f       <- y.f[[2]]
   num.pos <- sum(y>0)

   if (length(f) <= 1 || num.pos == 0 || num.pos == length(f))
   {
      return (-1.0)
   }
   else
   {
      return (gbm.roc.area(obs=y, pred=f))
   }
}

ir.measure.mrr <- function(y.f, max.rank)
{
   y       <- y.f[[1]]
   f       <- y.f[[2]]
   num.pos <- sum(y>0)

   if (length(f) <= 1 || num.pos == 0 || num.pos == length(f))
   {
      return (-1.0)
   }

   ord         <- order(f, decreasing=TRUE)
   min.idx.pos <- min(which(y[ord]>0))

   if (min.idx.pos <= max.rank)
   {
      return (1.0 / min.idx.pos)
   }
   else
   {
      return (0.0)
   }
}

ir.measure.map <- function(y.f, max.rank=0)
{
   # Note: max.rank is meaningless for MAP

   y         <- y.f[[1]]
   f         <- y.f[[2]]
   ord       <- order(f, decreasing=TRUE)
   idx.pos   <- which(y[ord]>0)
   num.pos   <- length(idx.pos)

   if (length(f) <= 1 || num.pos == 0 || num.pos == length(f))
   {
      return (-1.0)
   }

   # Above and including the rank of the i-th positive result,
   # there are exactly i positives and rank(i) total results
   return (sum((1:length(idx.pos))/idx.pos) / num.pos)
}

ir.measure.ndcg <- function(y.f, max.rank)
{
   y         <- y.f[[1]]
   f         <- y.f[[2]]

   if (length(f) <= 1 || all(diff(y)==0))
   {
      return (-1.0)
   }

   num.items <- min(length(f), max.rank)
   ord       <- order(f, decreasing=TRUE)

   dcg       <- sum(y[ord][1:num.items] / log2(2:(num.items+1)))

   # The best possible DCG: order by target
   ord.max   <- order(y, decreasing=TRUE)
   dcg.max   <- sum(y[ord.max][1:num.items] / log2(2:(num.items+1)))

   # Normalize
   return (dcg / dcg.max)
}
