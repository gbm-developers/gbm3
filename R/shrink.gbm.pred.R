shrink.gbm.pred <- function(object,newdata,n.trees,
                            lambda=rep(1,length(object$var.names)),
                            ...)
{
   if(length(lambda) != length(object$var.names))
   {
      stop("lambda must have the same length as the number of variables in the gbm object.")
   }

   if(!is.null(object$Terms))
   {
      x <- model.frame(delete.response(object$Terms),
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
         j <- match(levels(x[,i]), object$var.levels[[i]])
         if(any(is.na(j)))
         {
            stop(paste("New levels for variable ",
                       object$var.names[i],": ",
                       levels(x[,i])[is.na(j)],sep=""))
         }
         x[,i] <- as.numeric(x[,i])-1
      }
   }

   x <- as.vector(unlist(x))
   if(missing(n.trees) || any(n.trees > object$n.trees))
   {
      n.trees <- n.trees[n.trees<=object$n.trees]
      if(length(n.trees)==0) n.trees <- object$n.trees
      warning("n.trees not specified or some values exceeded number fit so far. Using ",n.trees,".")
   }
   # sort n.trees so that predictions are easier to generate and store
   n.trees <- sort(n.trees)

   predF <- .Call("gbm_shrink_pred",
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

   return(predF)
}
