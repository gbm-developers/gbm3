gbm.loss <- function(y, f, w, offset, dist, baseline, group=NULL, max.rank=NULL)
{
   if (!is.na(offset))
   {
      f <- offset+f
   }

   if (dist$name != "pairwise")
   {
      switch(dist$name,
             gaussian = weighted.mean((y - f)^2,w) - baseline,
             bernoulli = -2*weighted.mean(y*f - log(1+exp(f)),w) - baseline,
             laplace = weighted.mean(abs(y-f),w) - baseline,
             adaboost = weighted.mean(exp(-(2*y-1)*f),w) - baseline,
             poisson = -2*weighted.mean(y*f-exp(f),w) - baseline,
             stop(paste("Distribution",dist$name,"is not yet supported for method=permutation.test.gbm")))
   }
   else # dist$name == "pairwise"
   {
      if (is.null(dist$metric))
      {
         stop("No metric specified for distribution 'pairwise'")
      }
      if (!is.element(dist$metric, c("conc", "ndcg", "map", "mrr")))
      {
         stop("Invalid metric '", dist$metric, "' specified for distribution 'pairwise'")
      }
      if (is.null(group))
      {
         stop("For distribution 'pairwise', parameter 'group' has to be supplied")
      }
      # Loss = 1 - utility
      (1 - perf.pairwise(y, f, group, dist$metric, w, max.rank)) - baseline
   }
}
