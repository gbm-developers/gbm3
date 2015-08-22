#' Predict method for GBM Model Fits
#' 
#' Predicted values based on a generalized boosted model object
#' 
#' \code{predict.gbm} produces predicted values for each observation in
#' \code{newdata} using the first \code{n.trees} iterations of the boosting
#' sequence. If \code{n.trees} is a vector than the result is a matrix with
#' each column representing the predictions from gbm models with
#' \code{n.trees[1]} iterations, \code{n.trees[2]} iterations, and so on.
#' 
#' The predictions from \code{gbm} do not include the offset term. The user may
#' add the value of the offset to the predicted value if desired.
#' 
#' If \code{object} was fit using \code{\link{gbm.fit}}, there will be no
#' \code{Terms} component. Therefore, the user has greater responsibility to
#' make sure that \code{newdata} is of the same format (order and number of
#' variables) as the one originally used to fit the model.
#' 
#' @param object Object of class inheriting from (\code{\link{gbm.object}})
#' @param newdata Data frame of observations for which to make predictions
#' @param n.trees Number of trees used in the prediction. If \code{n.trees} is
#' a vector, predictions are returned for each iteration specified.
#' @param type The scale on which gbm makes the predictions
#' @param single.tree If \code{single.tree=TRUE} then \code{predict.gbm}
#' returns only the predictions from tree(s) \code{n.trees}
#' @param \dots further arguments passed to or from other methods
#' @return Returns a vector of predictions. By default the predictions are on
#' the scale of f(x). For example, for the Bernoulli loss the returned value is
#' on the log odds scale, poisson loss on the log scale, and coxph is on the
#' log hazard scale.
#' 
#' If \code{type="response"} then \code{gbm} converts back to the same scale as
#' the outcome. Currently the only effect this will have is returning
#' probabilities for bernoulli and expected counts for poisson. For the other
#' distributions "response" and "link" return the same.
#' @author Greg Ridgeway \email{gregridgeway@@gmail.com}
#' @seealso \code{\link{gbm}}, \code{\link{gbm.object}}
#' @keywords models regression
#' @export
predict.gbm <- function(object,newdata,n.trees,
                        type="link",
                        single.tree = FALSE,
                        ...)
{
   if ( missing( newdata ) ){
      newdata <- reconstructGBMdata(object)
   }
   if (missing(n.trees)){
      if ( object$train.fraction < 1 ){
         n.trees <- gbm.perf( object, method="test", plot.it = FALSE )
      } else if (!is.null(object$cv.error)){
         n.trees <- gbm.perf( object, method="cv", plot.it = FALSE )
      } else{
        n.trees <- length(object$train.error)
      }
      message(paste("Using", n.trees, "trees...\n"))
   } else if (length(n.trees) == 0){
      stop("n.trees cannot be NULL or a vector of zero length")
   }

   if(!is.element(type, c("link","response" )))
   {
      stop("type must be either 'link' or 'response'")
   }
   if(!is.null(object$Terms))
   {
      x <- model.frame(terms(reformulate(object$var.names)),
                       newdata,
                       na.action=na.pass)
   }
   else
   {
      x <- newdata
   }

   cRows <- nrow(x)
   cCols <- ncol(x)

   for(i in 1:cCols)
   {
      if(is.factor(x[,i]))
      {
        if (length(levels(x[,i])) > length(object$var.levels[[i]])) {
          new.compare <- levels(x[,i])[1:length(object$var.levels[[i]])]
        } else {
          new.compare <- levels(x[,i])
        }
        if (!identical(object$var.levels[[i]], new.compare)) {
          x[,i] <- factor(x[,i], union(object$var.levels[[i]], levels(x[,i])))
        }
        x[,i] <- as.numeric(x[,i])-1
      }
   }

   x <- as.vector(unlist(x, use.names=FALSE))
   if(missing(n.trees) || any(n.trees > object$n.trees))
   {
      n.trees[n.trees>object$n.trees] <- object$n.trees
      warning("Number of trees not specified or exceeded number fit so far. Using ",paste(n.trees,collapse=" "),".")
   }
   i.ntree.order <- order(n.trees)

   # Next if block for compatibility with objects created with version 1.6.
   if (is.null(object$num.classes)){
       object$num.classes <- 1
   }

   predF <- .Call("gbm_pred",
                  X=matrix(x, cRows, cCols),
                  n.trees=as.integer(n.trees[i.ntree.order]),
                  initF=object$initF,
                  trees=object$trees,
                  c.split=object$c.split,
                  var.type=as.integer(object$var.type),
                  single.tree = as.integer(single.tree),
                  PACKAGE = "gbm")

   if((length(n.trees) > 1) || (object$num.classes > 1))
   {
      predF <- matrix(predF, ncol=length(n.trees), byrow=FALSE)
      colnames(predF) <- n.trees
      predF[,i.ntree.order] <- predF
   }

   if(type=="response")
   {
      if(is.element(object$distribution$name, c("bernoulli", "pairwise")))
      {
         predF <- 1/(1+exp(-predF))
      } else
      if(object$distribution$name=="poisson")
      {
         predF <- exp(predF)
      } else
      if(object$distribution$name=="gamma")
      {
         predF <- exp(predF)
      } else
      if(object$distribution$name=="tweedie")
      {
         predF <- exp(predF)
      }
      else if (object$distribution$name == "adaboost"){
         predF <- 1 / (1 + exp(-2*predF))
      }
      if(length(n.trees)==1)
      {
         predF <- as.vector(predF)
      }
   }

   if(!is.null(attr(object$Terms,"offset")))
   {
      warning("predict.gbm does not add the offset to the predicted values.")
   }

   return(predF)
}
