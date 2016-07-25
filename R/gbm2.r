#' GBM2
#' 
#' Fits generalized boosted regression models - new API. This prepares the inputs, performing tasks
#' such as creating cv folds, before calling \code{gbm2.fit} to call the underlying C++ and fit a generalized
#' boosting model.
#' 
#' @usage  gbm2(formula, distribution=gbm_dist("Gaussian", ...), data, weights=rep(1, nrow(data)), offset=rep(0, nrow(data)),
#' train_params=training_params(num_trees=100, interaction_depth=1, min_num_obs_in_node=10, 
#' shrinkage=0.001, bag_fraction=0.5, id=seq_len(nrow(data)), num_train=round(0.5 * nrow(data)), num_features=ncol(data)-1) ,
#' var_monotone=NULL, var_names=NULL, cv_folds=1, cv_class_stratify=FALSE, fold_id=NULL, keep_gbm_data=FALSE, is_verbose=FALSE)
#' 
#' @param formula a symbolic description of the model to be fit.  The formula may include
#' an offset term (e.g. y~offset(n) + x).
#' 
#' @param distribution a GBMDist object specifying the distribution and any additional parameters needed.
#' 
#' @param data a data frame containing the variables in the model.  By default, the variables are taken from the 
#' environment. 
#' 
#' @param weights optional vector of weights used in the fitting process.  These weights must be positive but 
#' need not be normalized. By default they are set to 1 for each data row.
#' 
#' @param offset  optional vector specifying the model offset; must be positive.  This defaults to a vector of 0's, the length
#' of which is equal to the number rows of data.
#' 
#' @param train_params  a GBMTrainParams object which specifies the parameters used in growing decision trees.
#' 
#' @param var_monotone optional vector, the same length as the number of predictors, indicating the relationship
#' each variable has with the outcome.  It have a monotone increasing (+1) or decreasing (-1) or an arbitrary relationship.
#' 
#' @param var_names a vector of strings of containing the names of the predictor variables.
#' 
#' @param keep_gbm_data a bool specifying whether or not the gbm_data object created in this method should be
#' stored in the results.
#' 
#' @param cv_folds a positive integer specifying the number of folds to be used in cross-validation of the gbm
#' fit. If cv_folds > 1 then cross-validation is performed; the default of cv_folds is 1.
#' 
#' @param cv_class_stratify a bool specifying whether or not to stratify via response outcome. Currently only 
#' applies to "Bernoulli" distribution and defaults to false.
#' 
#' @param fold_id An optional vector of values identifying what fold
#' each observation is in. If supplied, cv_folds can be missing. Note:
#' Multiple rows of the same observation must have the same fold_id.
#' 
#' @param is_verbose if TRUE, gbm2 will print out progress and performance of the fit.
#' 
#' @return a \code{GBMFit} object.
#' 
#' @export gbm2
#' 

gbm2 <- function(formula, distribution=gbm_dist("Gaussian", ...), data, weights=rep(1, nrow(data)), offset=rep(0, nrow(data)),
                 train_params=training_params(num_trees=100, interaction_depth=1, min_num_obs_in_node=10, 
                 shrinkage=0.001, bag_fraction=0.5, id=seq(nrow(data)), num_train=nrow(data), num_features=ncol(data)-1), 
                 var_monotone=NULL, var_names=NULL,  cv_folds=1, cv_class_stratify=FALSE, fold_id=NULL,
                 keep_gbm_data=FALSE, is_verbose=FALSE) {
  
  # Extract the model
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "weights", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.pass
  mf[[1]] <- as.name("model.frame")
  m <- mf
  mf <- eval(mf, parent.frame())
  Terms <- attr(mf, "terms")
  y <- model.response(mf)
  w <- model.weights(mf)
  offset_mf <- model.offset(mf)
  
  # Set offset/ weights off defaults if specified
  if(!is.null(w))
    weights <- w
  if(!is.null(offset_mf))
    offset <- offset_mf
  
  # get the character name of the response variable
  var_names <- attributes(Terms)$term.labels
  x <- model.frame(terms(reformulate(var_names)),
                   data,
                   na.action=na.pass)
  
  # Check and infer folds if necessary
  check_cv_parameters(cv_folds, cv_class_stratify, fold_id, train_params)
  if (!is.null(fold_id)) {
    if (length(fold_id) != nrow(x)){
      stop("length(fold_id) inequal to number of rows.")
    }
    num_inferred_folds <- length(unique(fold_id))
    if (cv_folds != num_inferred_folds) {
      # Warn if cv_folds and fold_id disagree, but take fold_id.
      warning("CV folds changed from ", cv_folds, " to ", num_inferred_folds,
              " because of levels in fold_id.")
    } 
    cv_folds <- num_inferred_folds
    
    # Set fold_id from whatever it is to an integer ascending from 1. Lazy way.
    fold_id <- as.numeric(as.factor(fold_id))
  }
  
  # Call gbm2.fit 
  gbm_fit_obj <- gbm2.fit(x, y, distribution, weights, offset,
                         train_params, var_monotone, var_names, keep_gbm_data, cv_folds,
                         cv_class_stratify, fold_id, is_verbose)

  # Wrap up extra pieces 
  gbm_fit_obj$model <- m
  gbm_fit_obj$Terms <- Terms
 
  return(gbm_fit_obj)
}