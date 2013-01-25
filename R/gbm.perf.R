gbm.perf <- function(object,
                     plot.it=TRUE,
                     oobag.curve=FALSE,
                     overlay=TRUE,
                     method)
{
   smoother <- NULL

   if ( missing( method ) ){
      if ( object$train.fraction < 1 ){
         method <- "test"
      }
      else if ( !is.null( object$cv.error ) ){
         method <- "cv"
      }
      else { method <- "OOB" }
      cat( paste( "Using", method, "method...\n" ) )
   }

   if((method == "OOB") || oobag.curve)
   {
      if(object$bag.fraction==1)
      stop("Cannot compute OOB estimate or the OOB curve when bag.fraction=1")
      if(all(!is.finite(object$oobag.improve)))
      stop("Cannot compute OOB estimate or the OOB curve. No finite OOB estimates of improvement")
      x <- 1:object$n.trees
      smoother <- loess(object$oobag.improve~x,
                        enp.target=min(max(4,length(x)/10),50))
      smoother$y <- smoother$fitted
      smoother$x <- x

      best.iter.oob <- x[which.min(-cumsum(smoother$y))]
      best.iter <- best.iter.oob
   }

   if(method == "OOB")
   {
      warning("OOB generally underestimates the optimal number of iterations although predictive performance is reasonably competitive. Using cv.folds>0 when calling gbm usually results in improved predictive performance.")
   }

   if(method == "test")
   {
      best.iter.test <- which.min(object$valid.error)
      best.iter <- best.iter.test
   }

   if(method == "cv")
   {
      if(is.null(object$cv.error))
      stop("In order to use method=\"cv\" gbm must be called with cv.folds>1.")
      if(length(object$cv.error) < object$n.trees)
      warning("cross-validation error is not computed for any additional iterations run using gbm.more().")
      best.iter.cv <- which.min(object$cv.error)
      best.iter <- best.iter.cv
   }

   if(!is.element(method,c("OOB","test","cv")))
   stop("method must be cv, test, or OOB")

   if(plot.it)
   {
      par(mar=c(5,4,4,4)+.1)
      if (object$distribution$name !="pairwise")
      {
         ylab <- switch(substring(object$distribution$name,1,2),
                        ga="Squared error loss",
                        be="Bernoulli deviance",
                        po="Poisson deviance",
                        ad="AdaBoost exponential bound",
                        co="Cox partial deviance",
                        la="Absolute loss",
                        qu="Quantile loss",
                        mu="Multinomial deviance",
                        td="t-distribution deviance"
                        )
      }
      else # object$distribution$name =="pairwise"
      {
         ylab <- switch(object$distribution$metric,
                        conc ="Fraction of concordant pairs",
                        ndcg="Normalized discounted cumulative gain",
                        map ="Mean average precision",
                        mrr ="Mean reciprocal rank"
                        )
      }

      if(object$train.fraction==1)
      {  # HS Next line changed to scale axis to include other error
         #         ylim <- range(object$train.error)
         if ( method=="cv" ){ ylim <- range(object$train.error, object$cv.error) }
         else if ( method == "test" ){ ylim <- range( object$train.error, object$valid.error) }
         else { ylim <- range(object$train.error) }
      }
      else
      {
         ylim <- range(object$train.error,object$valid.error)
      }

      plot(object$train.error,
           ylim=ylim,
           type="l",
           xlab="Iteration",ylab=ylab)

      if(object$train.fraction!=1)
      {
         lines(object$valid.error,col="red")
      }
      if(method=="cv")
      {
         lines(object$cv.error,col="green")
      }
      if(!is.na(best.iter)) abline(v=best.iter,col="blue",lwd=2,lty=2)
      if(oobag.curve)
      {
         if(overlay)
         {
            par(new=TRUE)
            plot(smoother$x,
                 cumsum(smoother$y),
                 col="blue",
                 type="l",
                 xlab="",ylab="",
                 axes=FALSE)
            axis(4,srt=0)
            at <- mean(range(smoother$y))
            mtext(paste("OOB improvement in",ylab),side=4,srt=270,line=2)
            abline(h=0,col="blue",lwd=2)
         }

         plot(object$oobag.improve,type="l",
              xlab="Iteration",
              ylab=paste("OOB change in",ylab))
         lines(smoother,col="red",lwd=2)
         abline(h=0,col="blue",lwd=1)

         abline(v=best.iter,col="blue",lwd=1)
      }
   }

   return(best.iter)
}


perf.pairwise <- function(y, f, group, metric="ndcg", w=NULL, max.rank=0)
{
   func.name <- switch(metric,
                       conc = "ir.measure.conc",
                       mrr  = "ir.measure.mrr",
                       map  = "ir.measure.map",
                       ndcg = "ir.measure.ndcg",
                       stop(paste("Metric",metric,"is not supported"))
                       )

   # Optimization: for binary targets,
   # AUC is equivalent but faster than CONC
   if (metric == "conc" && all(is.element(y, 0:1)))
   {
      func.name <- "ir.measure.auc"
   }

   # Max rank = 0 means no cut off
   if (max.rank <= 0)
   {
      max.rank <- length(y)+1
   }

   # Random tie breaking in case of duplicate scores.
   # (Without tie breaking, we would overestimate if instances are
   # sorted descending on target)
   f <- f + 1E-10 * runif(length(f), min=-0.5, max=0.5)

   measure.by.group <- as.matrix(by(list(y, f), INDICES=group, FUN=get(func.name), max.rank=max.rank))

   # Exclude groups with single result or only negative or positive instances
   idx <- which((!is.null(measure.by.group)) & measure.by.group >= 0)

   if (is.null(w))
   {
      return (mean(measure.by.group[idx]))
   }
   else
   {
      # Assumption: weights are constant per group
      w.by.group <- tapply(w, group, mean)
      return (weighted.mean(measure.by.group[idx], w=w.by.group[idx]))
   }
}
