#' Methods for estimating relative influence
#' 
#' Helper functions for computing the relative influence of each variable in
#' the gbm object.
#' 
#' These functions offer the different
#' methods for computing the relative influence in \code{\link{summary.gbm}}.
#'
#' @aliases relative.influence permutation.test.gbm
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
#' @return By default, returns an unprocessed vector of estimated relative
#' influences. If the \code{scale.} and \code{sort.} arguments are used,
#' returns a processed version of the same.
#' @author Greg Ridgeway \email{gregridgeway@@gmail.com}
#' @details \code{\link{relative.influence}} is the same as that
#' described in Friedman (2001).
#' \code{\link{permutation.test.gbm}} randomly permutes each
#' predictor variable at a time and computes the associated reduction in
#' predictive performance. This is similar to the variable importance measures
#' Breiman uses for random forests, but \code{gbm} currently computes using the
#' entire training dataset (not the out-of-bag observations).

#' @seealso \code{\link{summary.gbm}}
#' @references J.H. Friedman (2001). "Greedy Function Approximation: A Gradient
#' Boosting Machine," Annals of Statistics 29(5):1189-1232.
#' 
#' L. Breiman (2001).
#' \href{http://oz.berkeley.edu/users/breiman/randomforest2001.pdf}{Random
#' Forests}.
#' @keywords hplot
#' @export
relative.influence <- function(object, n.trees, scale. = FALSE, sort. = FALSE){
   if( missing( n.trees ) ){
      if ( object$train.fraction < 1 ){
         n.trees <- gbm.perf( object, method="test", plot.it=FALSE )
      }
      else if ( !is.null( object$cv.error ) ){
         n.trees <- gbm.perf( object, method="cv", plot.it = FALSE )
      }
      else{
         n.trees <- object$n.trees
      }
      message(paste( "n.trees not given. Using", n.trees, "trees.\n" ))
      
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

   if (scale.) rel.inf <- rel.inf / max(rel.inf)
   if (sort.) rel.inf <- rev(sort(rel.inf))

   rel.inf
}
