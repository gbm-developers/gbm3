pretty.gbm.tree <- function(object,i.tree=1)
{
   if((i.tree<1) || (i.tree>length(object$trees)))
   {
      stop("i.tree is out of range. Must be less than ",length(object$trees))
   }
   else
   {
      temp <- data.frame(object$trees[[i.tree]])
      names(temp) <- c("SplitVar","SplitCodePred","LeftNode",
                       "RightNode","MissingNode","ErrorReduction",
                       "Weight","Prediction")
      row.names(temp) <- 0:(nrow(temp)-1)
   }
   return(temp)
}
