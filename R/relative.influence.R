relative.influence <- function(object,
                               n.trees,
                               scale. = FALSE,
                               sort. = FALSE )
{

   if( missing( n.trees ) ){
      if ( object$train.fraction < 1 ){
         n.trees <- gbm.perf( object, method="test", plot.it=FALSE )
      }
      else if ( !is.null( object$cv.error ) ){
         n.trees <- gbm.perf( object, method="cv", plot.it = FALSE )
      }
      else{
         # If dist=multinomial, object$n.trees = n.trees * num.classes
         # so use the following instead.
         n.trees <- length( object$train.error )
      }
      cat( paste( "n.trees not given. Using", n.trees, "trees.\n" ) )
      if (object$distribution == "multinomial"){
          n.trees <- n.trees * object$num.classes
      }
   }
   get.rel.inf <- function(obj)
   {
      lapply(split(obj[[6]],obj[[1]]),sum) # 6 - Improvement, 1 - var name
   }

   temp <- unlist(lapply(object$trees[1:n.trees],get.rel.inf))
   rel.inf.compact <- unlist(lapply(split(temp,names(temp)),sum))
   rel.inf.compact <- rel.inf.compact[names(rel.inf.compact)!="-1"]

   # rel.inf.compact excludes those variable that never entered the model
   # insert 0's for the excluded variables
   rel.inf <- rep(0,length(object$var.names))
   i <- as.numeric(names(rel.inf.compact))+1
   rel.inf[i] <- rel.inf.compact

   names(rel.inf) <- object$var.names

   if (scale.){
      rel.inf <- rel.inf / max(rel.inf)
   }
   if (sort.){
      rel.inf <- rev(sort(rel.inf))
   }

   return(rel.inf=rel.inf)
}
