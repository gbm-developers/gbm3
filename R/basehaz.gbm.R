# compute Breslow estimator of the baseline hazard function
basehaz.gbm <- function(t,delta,f.x,
                        t.eval=NULL,
                        smooth=FALSE,
                        cumulative=TRUE)
{
   t.unique <- sort(unique(t[delta==1]))
   alpha <- length(t.unique)
   for(i in 1:length(t.unique))
   {
      alpha[i] <- sum(t[delta==1]==t.unique[i])/
                     sum(exp(f.x[t>=t.unique[i]]))
   }

   if(!smooth && !cumulative)
   {
      if(!is.null(t.eval))
      {
         stop("Cannot evaluate unsmoothed baseline hazard at t.eval.")
      }
   } else
   if(smooth && !cumulative)
   {
      lambda.smooth <- supsmu(t.unique,alpha)
   } else
   if(smooth && cumulative)
   {
      lambda.smooth <- supsmu(t.unique,cumsum(alpha))
   } else # (!smooth && cumulative) - THE DEFAULT
   {
      lambda.smooth <- list(x=t.unique,y=cumsum(alpha))
   }

   if(!is.null(t.eval))
   {
      obj <- approx(lambda.smooth$x,lambda.smooth$y,xout=t.eval)$y
   } else
   {
      obj <- approx(lambda.smooth$x,lambda.smooth$y,xout=t)$y
   }

   return(obj)
}
