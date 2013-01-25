permutation.test.gbm <- function(object,
                                 n.trees)
{
   # get variables used in the model
   i.vars <- sort(unique(unlist(lapply(object$trees[1:n.trees],
                                       function(x){unique(x[[1]])}))))
   i.vars <- i.vars[i.vars!=-1] + 1
   rel.inf <- rep(0,length(object$var.names))

   if(!is.null(object$data))
   {
      y            <- object$data$y
      os           <- object$data$offset
      Misc         <- object$data$Misc
      w            <- object$data$w
      x            <- matrix(object$data$x, ncol=length(object$var.names))
      object$Terms <- NULL # this makes predict.gbm take x as it is

      if (object$distribution$name == "pairwise")
      {
         # group and cutoff are only relevant for distribution "pairwise"
         # in this case, the last element specifies the max rank
         # max rank = 0 means no cut off
         group     <- Misc[1:length(y)]
         max.rank  <- Misc[length(y)+1]
      }
   }
   else
   {
      stop("Model was fit with keep.data=FALSE. permutation.test.gbm has not been implemented for that case.")
   }

   # the index shuffler
   j <- sample(1:nrow(x))
   for(i in 1:length(i.vars))
   {
      x[ ,i.vars[i]]  <- x[j,i.vars[i]]

      new.pred <- predict.gbm(object,newdata=x,n.trees=n.trees)
      rel.inf[i.vars[i]] <- gbm.loss(y,new.pred,w,os,
                                     object$distribution,
                                     object$train.error[n.trees],
                                     group,
                                     max.rank)

      x[j,i.vars[i]] <- x[ ,i.vars[i]]
   }

   return(rel.inf=rel.inf)
}
