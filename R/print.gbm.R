# print and summary functions for gbm

#' Print model summary
#' 
#' Display basic information about a \code{gbm} object.
#' 
#' Prints some information about the model object. In particular, this method
#' prints the call to \code{gbm()}, the type of loss function that was used,
#' and the total number of iterations.
#' 
#' If cross-validation was performed, the 'best' number of trees as estimated
#' by cross-validation error is displayed. If a test set was used, the 'best'
#' number of trees as estimated by the test set error is displayed.
#' 
#' The number of available predictors, and the number of those having non-zero
#' influence on predictions is given (which might be interesting in data mining
#' applications).
#' 
#' If multinomial, bernoulli or adaboost was used, the confusion matrix and
#' prediction accuracy are printed (objects being allocated to the class with
#' highest probability for multinomial and bernoulli). These classifications
#' are performed using the cross-validation fitted values.
#' 
#' If the 'distribution' was specified as gaussian, laplace, quantile or
#' t-distribution, a summary of the residuals is displayed.  The residuals are
#' the cross-validation residuals. Also, a pseudo R-squared value is displayed.
#' For Gaussian response, this is 1 - sum(r*r) / sum(z*z) where z = y -
#' mean(y). For the other distributions, this is 1 - (median(abs(r)) /
#' mad(y))^2, following the suggestion of Rousseeuw and Leroy (equation 3.11).
#' Note that this definition of a robust R-squared is contentious.
#' 
#' @param x an object of class \code{gbm}.
#' @param \dots arguments passed to \code{print.default}.
#' @author Harry Southworth, Daniel Edwards
#' @seealso \code{\link{gbm}}
#' @references P. J. Rousseeuw and A. M. Leroy, Robust Regression and Outlier
#' Detection, Wiley, 1987 (2003).
#' @keywords models nonlinear survival nonparametric
#' @examples
#' 
#' data(iris)
#' iris.mod <- gbm(Species ~ ., distribution="multinomial", data=iris,
#'                  n.trees=2000, shrinkage=0.01, cv.folds=5,
#'                  verbose=FALSE, n.cores=1)
#' iris.mod
#' #data(lung)
#' #lung.mod <- gbm(Surv(time, status) ~ ., distribution="coxph", data=lung,
#' #                 n.trees=2000, shrinkage=0.01, cv.folds=5,verbose =FALSE)
#' #lung.mod
#' @export
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

   if (is.null(x$cv.fitted)) {
       return(invisible())
   }
   
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

   invisible()
}



#' Summary of a gbm object
#' 
#' Computes the relative influence of each variable in the gbm object.
#' 
#' For \code{distribution="gaussian"} this returns exactly the reduction of
#' squared error attributable to each variable. For other loss functions this
#' returns the reduction attributable to each variable in sum of squared error
#' in predicting the gradient on each iteration. It describes the relative
#' influence of each variable in reducing the loss function. See the references
#' below for exact details on the computation.
#' 
#' @param object a \code{gbm} object created from an initial call to
#' \code{\link{gbm}}.
#' @param cBars the number of bars to plot. If \code{order=TRUE} the only the
#' variables with the \code{cBars} largest relative influence will appear in
#' the barplot. If \code{order=FALSE} then the first \code{cBars} variables
#' will appear in the plot. In either case, the function will return the
#' relative influence of all of the variables.
#' @param n.trees the number of trees used to generate the plot. Only the first
#' \code{n.trees} trees will be used.
#' @param plotit an indicator as to whether the plot is generated.
#' @param order an indicator as to whether the plotted and/or returned relative
#' influences are sorted.
#' @param method The function used to compute the relative influence.
#' \code{\link{relative.influence}} is the default and is the same as that
#' described in Friedman (2001). The other current (and experimental) choice is
#' \code{\link{permutation.test.gbm}}. This method randomly permutes each
#' predictor variable at a time and computes the associated reduction in
#' predictive performance. This is similar to the variable importance measures
#' Breiman uses for random forests, but \code{gbm} currently computes using the
#' entire training dataset (not the out-of-bag observations).
#' @param normalize if \code{FALSE} then \code{summary.gbm} returns the
#' unnormalized influence.
#' @param ...  other arguments passed to the plot function.
#' @return Returns a data frame where the first component is the variable name
#' and the second is the computed relative influence, normalized to sum to 100.
#' @author Greg Ridgeway \email{gregridgeway@@gmail.com}
#' @seealso \code{\link{gbm}}
#' @references J.H. Friedman (2001). "Greedy Function Approximation: A Gradient
#' Boosting Machine," Annals of Statistics 29(5):1189-1232.
#' 
#' L. Breiman (2001). \href{http://oz.berkeley.edu/users/breiman/randomforest2001.pdf}{Random Forests}.
#' @keywords hplot
#' @export
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
              xlab="Relative influence",
              las=1,...)
   }
   return(data.frame(var=object$var.names[i],
                     rel.inf=rel.inf[i]))
}
