# evaluates the objective function and gradient with respect to beta
# beta = log(lambda/(1-lambda))
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
