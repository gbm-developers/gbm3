#' Predict method for GBM Model Fits
#' 
#' Predicted values based on a generalized boosted model object - from
#' gbmt
#' 
#' \code{predict.GBMFit} produces predicted values for each
#' observation in a new dataset \code{newdata} using the first
#' \code{num_trees} iterations of the boosting sequence. If
#' \code{num_trees} is a vector than the result is a matrix with each
#' column representing the predictions from gbm models with
#' \code{num_trees[1]} iterations, \code{num_trees[2]} iterations, and
#' so on.
#' 
#' The predictions from \code{gbmt} do not include the offset
#' term. The user may add the value of the offset to the predicted
#' value if desired.
#' 
#' If \code{gbm_fit_obj} was fit using \code{\link{gbmt}}, there will
#' be no \code{Terms} component. Therefore, the user has greater
#' responsibility to make sure that \code{newdata} is of the same
#' format (order and number of variables) as the one originally used
#' to fit the model.
#' 
#' @param object Object of class inheriting from \code{GBMFit}.
#' @param newdata Data frame of observations for which to make
#' predictions
#' @param n.trees Number of trees used in the prediction. If
#' \code{n.trees} is a vector, predictions are returned for each
#' iteration specified.
#' @param type The scale on which gbm makes the predictions
#' @param single.tree If \code{single.tree=TRUE} then \code{gbm_predict}
#' returns only the predictions from tree(s) \code{n.trees}
#' @param \dots further arguments passed to or from other methods
#' @return Returns a vector of predictions. By default the predictions
#' are on the scale of f(x). For example, for the Bernoulli loss the
#' returned value is on the log odds scale, poisson loss on the log
#' scale, and coxph is on the log hazard scale.
#' 
#' If \code{type="response"} then \code{gbmt} converts back to the same scale as
#' the outcome. Currently the only effect this will have is returning
#' probabilities for bernoulli and expected counts for poisson. For the other
#' distributions "response" and "link" return the same.
#' @seealso \code{\link{gbmt}}
#' @keywords models regression
#' @importFrom stats predict
#' @export 
predict.GBMFit <- function(object, newdata, n.trees,
                           type="link", single.tree=FALSE,
                           ...)
{
  # Check inputs
  if(!is.element(type, c("link","response" ))) {
    stop("type must be either 'link' or 'response'")
  }
  
  if(missing(newdata) || !is.data.frame(newdata)) {
    stop("newdata must be provided as a data frame")
  }
  
  if(missing(n.trees)) {
    stop("Number of trees to be used in prediction must be provided.")
  }
  
  if (length(n.trees) == 0) {
    stop("n.trees cannot be NULL or a vector of zero length")
  }
  
  if(any(n.trees != as.integer(n.trees)) || is.na(all(n.trees == as.integer(n.trees)))
     || any(n.trees < 0)) {
    stop("n.trees must be a vector of non-negative integers")
  }
  
  if(!is.null(attr(object$Terms,"offset")))
  {
    warning("predict.GBMFit does not add the offset to the predicted values.")
  }
  
  # Get data
  if(!is.null(object$Terms)) {
    x <- model.frame(terms(reformulate(object$variables$var_names)),
                     newdata,
                     na.action=na.pass)
  } else {
    x <- newdata
  }
  
  # Convert predictor factors into appropriate numeric
  for(i in seq_len(ncol(x))) {
    if(is.factor(x[,i])) {
      if (length(levels(x[,i])) > length(object$variables$var_levels[[i]])) {
        new_compare <- levels(x[,i])[seq_len(length(object$variables$var_levels[[i]]))]
      } else {
        new_compare <- levels(x[,i])
      }
      
      if (!identical(object$variables$var_levels[[i]], new_compare)) {
        x[,i] <- factor(x[,i], union(object$variables$var_levels[[i]], levels(x[,i])))
      }
      
      x[,i] <- as.numeric(x[,i])-1
    }
  }
  
  if(any(n.trees > object$params$num_trees)) {
    n.trees[n.trees > object$params$num_trees] <- object$params$num_trees
    warning("Number of trees exceeded number fit so far. Using ", paste(n.trees,collapse=" "),".")
  }
  
  i.ntree.order <- order(n.trees)
  
  # Next if block for compatibility with objects created with 1.6
  if(is.null(object$num.classes)) {
    object$num.classes <- 1
  }
  
  predF <- .Call("gbm_pred",
                 X=as.matrix(as.data.frame(x)),
                 n.trees=as.integer(n.trees[order(n.trees)]),
                 initF=object$initF,
                 trees=trees(object),
                 c.split=object$c.split,
                 var.type=as.integer(object$variables$var_type),
                 single.tree = as.integer(single.tree),
                 PACKAGE = "gbm")
  
  # Convert into matrix of predictions
  if((length(n.trees) > 1) || (!is.null(object$num.classes) && (object$num.classes > 1))) {
    predF <- matrix(predF, ncol=length(n.trees), byrow=FALSE)
    colnames(predF) <- n.trees
    predF[, order(n.trees)] <- predF
  }
  
  # Adjust scale of predictions
  if(type=="response") {
    predF <- adjust_pred_scale(predF, object$distribution)
    if(length(n.trees)==1) {
      predF <- as.vector(predF)
    }
  }
  return(predF)
}
  
