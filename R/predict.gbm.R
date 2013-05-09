predict.gbm <- function(object,newdata,n.trees,
                        type="link",
                        single.tree = FALSE,
                        ...)
{
   if ( missing( newdata ) ){
      newdata <- reconstructGBMdata(object)
   }
   if ( missing(n.trees) ) {
      if ( object$train.fraction < 1 ){
         n.trees <- gbm.perf( object, method="test", plot.it = FALSE )
      }
      else if (!is.null(object$cv.error)){
         n.trees <- gbm.perf( object, method="cv", plot.it = FALSE )
      }
      else{ best <- length( object$train.error ) }
      cat( paste( "Using", n.trees, "trees...\n" ) )
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
                  X=as.double(x),
                  cRows=as.integer(cRows),
                  cCols=as.integer(cCols),
                  cNumClasses = as.integer(object$num.classes),
                  n.trees=as.integer(n.trees[i.ntree.order]),
                  initF=object$initF,
                  trees=object$trees,
                  c.split=object$c.split,
                  var.type=as.integer(object$var.type),
                  single.tree = as.integer(single.tree),
                  PACKAGE = "gbm")

   if((length(n.trees) > 1) || (object$num.classes > 1))
   {
      if(object$distribution$name=="multinomial")
      {
         predF <- array(predF, dim=c(cRows,object$num.classes,length(n.trees)))
         dimnames(predF) <- list(NULL, object$classes, n.trees)
         predF[,,i.ntree.order] <- predF
      } else
      {
         predF <- matrix(predF, ncol=length(n.trees), byrow=FALSE)
         colnames(predF) <- n.trees
         predF[,i.ntree.order] <- predF
      }
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
      }
      else if (object$distribution$name == "adaboost"){
         predF <- 1 / (1 + exp(-2*predF))
      }
      if(object$distribution$name=="multinomial")
      {
         pexp <- exp(predF)
         psum  <- apply(pexp,  c(1, 3), function(x) { x / sum(x) })
         # Transpose each 2d array
         predF <- aperm(psum, c(2, 1, 3))
      }

      if((length(n.trees)==1) && (object$distribution$name!="multinomial"))
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
