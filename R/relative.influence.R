#' Methods for estimating relative influence
#' 
#' Helper functions for computing the relative influence of each variable in
#' the gbm object.
#' 
#' This is not intended for end-user use. These functions offer the different
#' methods for computing the relative influence in \code{\link{summary.gbm}}.
#' \code{gbm.loss} is a helper function for \code{permutation.test.gbm}.
#' 
#' @aliases relative.influence permutation.test.gbm gbm.loss
#' @param object a \code{gbm} object created from an initial call to
#' \code{\link{gbm}}.
#' @param n.trees the number of trees to use for computations. If not provided,
#' the function will guess: if a test set was used in fitting, the number of
#' trees resulting in lowest test set error will be used; otherwise, if
#' cross-validation was performed, the number of trees resulting in lowest
#' cross-validation error will be used; otherwise, all trees will be used.
#' @param scale.  whether or not the result should be scaled. Defaults to
#' \code{FALSE}.
#' @param sort.  whether or not the results should be (reverse) sorted.
#' Defaults to \code{FALSE}.
#' @param y,f,w,offset,dist,baseline For \code{gbm.loss}: These components are
#' the outcome, predicted value, observation weight, offset, distribution, and
#' comparison loss function, respectively.
#' @param group,max.rank Used internally when \code{distribution =
#' \'pairwise\'}.
#' @return By default, returns an unprocessed vector of estimated relative
#' influences. If the \code{scale.} and \code{sort.} arguments are used,
#' returns a processed version of the same.
#' @author Greg Ridgeway \email{gregridgeway@@gmail.com}
#' @seealso \code{\link{summary.gbm}}
#' @references J.H. Friedman (2001). "Greedy Function Approximation: A Gradient
#' Boosting Machine," Annals of Statistics 29(5):1189-1232.
#' 
#' L. Breiman (2001).
#' \href{http://oz.berkeley.edu/users/breiman/randomforest2001.pdf}{Random
#' Forests}.
#' @keywords hplot
#' @export
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
      message(paste( "n.trees not given. Using", n.trees, "trees.\n" ))
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
