# print, show and summary functions for gbm

print.gbm <- function(x, ... ){
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

   #############################################################################

   d <- reconstructGBMdata(x)
   if (x$distribution$name == "multinomial"){
       n.class <- x$num.classes

       yn <- as.numeric(d[, x$response.name])
       
       p <- predict(x, n.trees=best , type = "response", newdata=d)
       p <- apply(p, 1, function(x , labels ){ labels[x == max(x)] }, labels = colnames(p))
       p <- as.numeric(as.factor(p))
       r <- yn

       conf.mat <- matrix(table(c(r + n.class * p, (n.class + (1:(n.class^2))))), 
                          nrow = n.class)
       conf.mat <- conf.mat - 1
       pred.acc <- round(100 * sum(diag(conf.mat)) / sum(conf.mat),2)
       conf.mat <- cbind(conf.mat, round(100*diag(conf.mat)/rowSums(conf.mat),2))
       dimnames(conf.mat) <- list(x$classes, c(x$classes, "Pred. Acc."))

       cat("\nConfusion matrix:\n")
       print(conf.mat)

       cat("\nPrediction Accuracy = ", pred.acc, "%\n", sep = "") 
   }
   else if (x$distribution$name %in% c("bernoulli", "adaboost", "huberized")){

       p <- predict( x , newdata=d, n.tree=best , type = "response")
       p <- ifelse( p < .5, 0, 1 )

       conf.mat <- matrix( table( c( d[,x$response.name] + 2 * p , 0:3 )), ncol=2 )
       conf.mat <- conf.mat - 1

       pred.acc <- round(100 * sum(diag(conf.mat)) / sum(conf.mat),2)

       conf.mat <- cbind(conf.mat,  round(100*diag(conf.mat)/rowSums(conf.mat),2))
       dimnames(conf.mat) <- list(c("0","1"), c("0", "1", "Pred. Acc."))

       cat("\nConfusion matrix:\n")
       print(conf.mat)

       cat("\nPrediction Accuracy = ", pred.acc, "%\n", sep = "")
   }
   else if ( x$distribution$name %in% c( "gaussian", "laplace", "poisson", "quantile", "bisquare", "tdist" ) ){
       r <- d[, 1] - predict( x, type="response", newdata=d, n.tree=best )
       if ( x$distribution$name == "poisson" ){
           cat( "Summary of response residuals:\n" )
       }
       else {
           cat( "Summary of residuals:\n" )
       }
       print( quantile( r ) )
       cat( "\n" )
   }
   

   #############################################################################
   
   invisible()
}

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
