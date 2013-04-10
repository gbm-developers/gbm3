# print, show and summary functions for gbm

print.gbm <- function(x, ... )
{
   if (!is.null(x$call)){ print(x$call) }
   dist.name <- x$distribution$name
   if (dist.name == "pairwise")
   {
      if (!is.null(x$distribution$max.rank) && x$distribution$max.rank > 0)
      {
           dist.name <- sprintf("pairwise (metric=%s, max.rank=%d)", x$distribution$metric,  x$distribution$max.rank)
      }
      else
      {
           dist.name <- sprintf("pairwise (metric=%s)", x$distribution$metric)
      }
   }
   cat( paste( "A gradient boosted model with", dist.name, "loss function.\n" ))
   cat( paste( length( x$train.error ), "iterations were performed.\n" ) )
   best <- length( x$train.error )
   if ( !is.null( x$cv.error ) )
   {
      best <- gbm.perf( x, plot.it = FALSE, method="cv" )
      cat( paste("The best cross-validation iteration was ", best, ".\n", sep = "" ) )
   }
   if ( x$train.fraction < 1 )
   {
      best <- gbm.perf( x, plot.it = FALSE, method="test" )
      cat( paste("The best test-set iteration was ", best, ".\n", sep = "" ) )
   }
   if ( is.null( best ) )
   {
      best <- length( x$train.error )
   }
   ri <- relative.influence( x, n.trees=best )
   cat( "There were", length( x$var.names ), "predictors of which",
       sum( ri > 0 ), "had non-zero influence.\n" )

   invisible()
}

show.gbm <- print.gbm

summary.gbm <- function(object,
                        cBars=length(object$var.names),
                        n.trees=object$n.trees,
                        plotit=TRUE,
                        order=TRUE,
                        method=relative.influence,
                        normalize=TRUE,
                        ...)
{
   if(n.trees < 1)
   {
      stop("n.trees must be greater than 0.")
   }
   if(n.trees > object$n.trees)
   {
      warning("Exceeded total number of GBM terms. Results use n.trees=",object$n.trees," terms.\n")
      n.trees <- object$n.trees
   }

   rel.inf <- method(object,n.trees)
   rel.inf[rel.inf<0] <- 0

   if(order)
   {
      i <- order(-rel.inf)
   }
   else
   {
      i <- 1:length(rel.inf)
   }
   if(cBars==0) cBars <- min(10,length(object$var.names))
   if(cBars>length(object$var.names)) cBars <- length(object$var.names)

   if(normalize) rel.inf <- 100*rel.inf/sum(rel.inf)

   if(plotit)
   {
      barplot(rel.inf[i[cBars:1]],
              horiz=TRUE,
              col=rainbow(cBars,start=3/6,end=4/6),
              names=object$var.names[i[cBars:1]],
              xlab="Relative influence",...)
   }
   return(data.frame(var=object$var.names[i],
                     rel.inf=rel.inf[i]))
}
