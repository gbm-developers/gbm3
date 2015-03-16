# evaluates the objective function and gradient with respect to beta
# beta = log(lambda/(1-lambda))


#' L1 shrinkage of the predictor variables in a GBM
#' 
#' Performs recursive shrinkage in each of the trees in a GBM fit using
#' different shrinkage parameters for each variable.
#' 
#' This function is currently experimental. Used in conjunction with a gradient
#' ascent search for inclusion of variables.
#' 
#' @param object A \code{\link{gbm.object}}
#' @param n.trees the number of trees to use
#' @param lambda a vector with length equal to the number of variables
#' containing the shrinkage parameter for each variable
#' @param \dots other parameters (ignored)
#' @return \item{predF}{Predicted values from the shrunken tree}
#' \item{objective}{The value of the loss function associated with the
#' predicted values} \item{gradient}{A vector with length equal to the number
#' of variables containing the derivative of the objective function with
#' respect to beta, the logit transform of the shrinkage parameter for each
#' variable}
#' @section Warning: This function is experimental.
#' @author Greg Ridgeway \email{gregridgeway@@gmail.com}
#' @seealso \code{\link{shrink.gbm.pred}}, \code{\link{gbm}}
#' @references Hastie, T. J., and Pregibon, D.
#' \href{http://www-stat.stanford.edu/~hastie/Papers/shrinktree.ps}{Shrinking
#' Trees}. AT&T Bell Laboratories Technical Report (March 1990).
#' @keywords methods
shrink.gbm <- function(object,n.trees,
                       lambda=rep(10,length(object$var.names)),
                       ...)
{
   if(length(lambda) != length(object$var.names))
   {
      stop("lambda must have the same length as the number of variables in the gbm object.")
   }

   if(is.null(object$data))
   {
      stop("shrink.gbm requires keep.data=TRUE when gbm model is fit.")
   }

   y <- object$data$y
   x <- object$data$x

   cCols <- length(object$var.names)
   cRows <- length(x)/cCols


   if(missing(n.trees) || (n.trees > object$n.trees))
   {
      n.trees <- object$n.trees
      warning("n.trees not specified or some values exceeded number fit so far. Using ",n.trees,".")
   }

   result <- .Call("gbm_shrink_gradient",
                   y=as.double(y),
                   X=as.double(x),
                   cRows=as.integer(cRows),
                   cCols=as.integer(cCols),
                   n.trees=as.integer(n.trees),
                   initF=object$initF,
                   trees=object$trees,
                   c.split=object$c.split,
                   var.type=as.integer(object$var.type),
                   depth=as.integer(object$interaction.depth),
                   lambda=as.double(lambda),
                   PACKAGE = "gbm")

   names(result) <- c("predF","objective","gradient")

   return(result)
}
