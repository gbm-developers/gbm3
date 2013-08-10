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

       p <- apply(x$cv.fitted, 1, function(x , labels){ labels[x == max(x)] }, labels = x$classes)
       p <- as.numeric(as.factor(p))
       r <- yn

       conf.mat <- matrix(table(c(r + n.class * p, (n.class + (1:(n.class^2))))),
                          nrow = n.class)
       conf.mat <- conf.mat - 1
       pred.acc <- round(100 * sum(diag(conf.mat)) / sum(conf.mat),2)
       conf.mat <- cbind(conf.mat, round(100*diag(conf.mat)/rowSums(conf.mat),2))
       dimnames(conf.mat) <- list(x$classes, c(x$classes, "Pred. Acc."))

       cat("\nCross-validation confusion matrix:\n")
       print(conf.mat)

       cat("\nCross-validation prediction Accuracy = ", pred.acc, "%\n", sep = "") 
   }
   else if (x$distribution$name %in% c("bernoulli", "adaboost", "huberized")){

       p <- 1 / (1 + exp(-x$cv.fitted))
       p <- ifelse(p < .5, 0, 1)

       conf.mat <- matrix(table(c(d[, x$response.name] + 2 * p , 0:3)), ncol=2)
       conf.mat <- conf.mat - 1

       pred.acc <- round(100 * sum(diag(conf.mat)) / sum(conf.mat),2)

       conf.mat <- cbind(conf.mat,  round(100*diag(conf.mat)/rowSums(conf.mat),2))
       dimnames(conf.mat) <- list(c("0","1"), c("0", "1", "Pred. Acc."))

       cat("\nCross-validation confusion matrix:\n")
       print(conf.mat)

       cat("\nCross-validation prediction Accuracy = ", pred.acc, "%\n", sep = "")
   }
   else if (x$distribution$name %in% c("gaussian", "laplace", "quantile", "tdist")){
       r <- d[, 1] - x$cv.fitted

       cat("\nSummary of cross-validation residuals:\n" )
       print(quantile(r))
       cat("\n")
       
       # Do pseudo R^2
       if (x$distribution$name == "gaussian"){
           yadj <- d[, 1] - mean(d[, 1])
           R2 <- 1 - sum(r^2)/sum(yadj^2)
           cat("Cross-validation pseudo R-squared: ", signif(R2, 3), "\n")
       }
       else { # Rousseeuw & Leroy, page 44
           R2 <- 1 - (median(abs(r)) / mad(d[, 1]))^2
           cat("Cross-validation robust pseudo R-squared: ", signif(R2, 3), "\n")
       }
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
