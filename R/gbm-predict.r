#' Predict method for GBM Model Fits 
#' 
#' Predicted values based on a generalized boosted model object - from gbmt
#' 
#' \code{gbm_predict} produces predicted values for each observation in a new_data_set
#' \code{new_data} using the first \code{num_trees} iterations of the boosting
#' sequence. If \code{num_trees} is a vector than the result is a matrix with
#' each column representing the predictions from gbm models with
#' \code{num_trees[1]} iterations, \code{num_trees[2]} iterations, and so on.
#' 
#' The predictions from \code{gbmt} do not include the offset term. The user may
#' add the value of the offset to the predicted value if desired.
#' 
#' If \code{gbm_fit_obj} was fit using \code{\link{gbmt}}, there will be no
#' \code{Terms} component. Therefore, the user has greater responsibility to
#' make sure that \code{new_data} is of the same format (order and number of
#' variables) as the one originally used to fit the model.
#' 
#' @param gbm_fit_obj Object of class inheriting from \code{GBMFit}.
#' @param new_data Data frame of observations for which to make predictions
#' @param num_trees Number of trees used in the prediction. If \code{num_trees} is
#' a vector, predictions are returned for each iteration specified.
#' @param type The scale on which gbm makes the predictions
#' @param single_tree If \code{single_tree=TRUE} then \code{gbm_predict}
#' returns only the predictions from tree(s) \code{num_trees}
#' @param \dots further arguments passed to or from other methods
#' @return Returns a vector of predictions. By default the predictions are on
#' the scale of f(x). For example, for the Bernoulli loss the returned value is
#' on the log odds scale, poisson loss on the log scale, and coxph is on the
#' log hazard scale.
#' 
#' If \code{type="response"} then \code{gbmt} converts back to the same scale as
#' the outcome. Currently the only effect this will have is returning
#' probabilities for bernoulli and expected counts for poisson. For the other
#' distributions "response" and "link" return the same.
#' @seealso \code{\link{gbmt}}
#' @keywords models regression
#' @export
#' 
predict <- function(gbm_fit_obj, new_data, num_trees,
                    type="link", single_tree=FALSE,
                    ...) {
  UseMethod("predict", gbm_fit_obj)
}

#' @name predict
#' @export 
predict.GBMFit <- function(gbm_fit_obj, new_data, num_trees,
                        type="link", single_tree=FALSE,
                        ...)
{
  # Check inputs
  if(!is.element(type, c("link","response" ))) {
    stop("type must be either 'link' or 'response'")
  }
  
  if(missing(new_data) || !is.data.frame(new_data)) {
    stop("new_data must be provided as a data frame")
  }
  
  if(missing(num_trees)) {
    stop("Number of trees to be used in prediction must be provided.")
  }
  
  if (length(num_trees) == 0) {
    stop("num_trees cannot be NULL or a vector of zero length")
  }
  
  if(any(num_trees != as.integer(num_trees)) || is.na(all(num_trees == as.integer(num_trees)))
     || any(num_trees < 1)) {
    stop("num_trees must be a vector of positive integers")
  }
  
  if(!is.null(attr(gbm_fit_obj$Terms,"offset")))
  {
    warning("predict.GBMFit does not add the offset to the predicted values.")
  }
  
  # Get data
  if(!is.null(gbm_fit_obj$Terms)) {
    x <- model.frame(terms(reformulate(gbm_fit_obj$variables$var_names)),
                     new_data,
                     na.action=na.pass)
  } else {
    x <- new_data
  }
  
  # Convert predictor factors into appropriate numeric
  for(i in seq_len(ncol(x))) {
    if(is.factor(x[,i])) {
      if (length(levels(x[,i])) > length(gbm_fit_obj$variables$var_levels[[i]])) {
        new_compare <- levels(x[,i])[seq_len(length(gbm_fit_obj$variables$var_levels[[i]]))]
      } else {
        new_compare <- levels(x[,i])
      }
      
      if (!identical(gbm_fit_obj$variables$var_levels[[i]], new_compare)) {
        x[,i] <- factor(x[,i], union(gbm_fit_obj$variables$var_levels[[i]], levels(x[,i])))
      }
      
      x[,i] <- as.numeric(x[,i])-1
    }
  }
  
  if(any(num_trees > gbm_fit_obj$params$num_trees)) {
    num_trees[num_trees > gbm_fit_obj$params$num_trees] <- gbm_fit_obj$params$num_trees
    warning("Number of trees exceeded number fit so far. Using ", paste(num_trees,collapse=" "),".")
  }
  
  i.ntree.order <- order(num_trees)
  
  # Next if block for compatibility with objects created with 1.6
  if(is.null(gbm_fit_obj$num.classes)) {
    gbm_fit_obj$num.classes <- 1
  }
  
  predF <- .Call("gbm_pred",
                 X=as.matrix(as.data.frame(x)),
                 n.trees=as.integer(num_trees[order(num_trees)]),
                 initF=gbm_fit_obj$initF,
                 trees=gbm_fit_obj$trees,
                 c.split=gbm_fit_obj$c.split,
                 var.type=as.integer(gbm_fit_obj$variables$var_type),
                 single.tree = as.integer(single_tree),
                 PACKAGE = "gbm")
  
  # Convert into matrix of predictions
  if((length(num_trees) > 1) || (!is.null(gbm_fit_obj$num.classes) && (gbm_fit_obj$num.classes > 1))) {
    predF <- matrix(predF, ncol=length(num_trees), byrow=FALSE)
    colnames(predF) <- num_trees
    predF[, order(num_trees)] <- predF
  }
  
  # Adjust scale of predictions
  if(type=="response") {
    predF <- adjust_pred_scale(predF, gbm_fit_obj$distribution)
    if(length(num_trees)==1) {
      predF <- as.vector(predF)
    }
  }
  return(predF)
}
  