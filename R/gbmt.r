##' GBMT
##' 
##' Fits generalized boosted regression models - new API. This
##' prepares the inputs, performing tasks such as creating cv folds,
##' before calling \code{gbmt_fit} to call the underlying C++ and fit
##' a generalized boosting model.
##' 
##' @param formula a symbolic description of the model to be fit.  The
##' formula may include an offset term (e.g. y~offset(n) + x).
##' 
##' @param distribution a \code{GBMDist} object specifying the
##' distribution and any additional parameters needed. If not
##' specified then the distribution will be guessed.
##' 
##' @param data a data frame containing the variables in the model.
##' By default, the variables are taken from the environment.
##' 
##' @param weights optional vector of weights used in the fitting
##' process.  These weights must be positive but need not be
##' normalized. By default they are set to 1 for each data row.
##' 
##' @param offset optional vector specifying the model offset; must be
##' positive.  This defaults to a vector of 0's, the length of which
##' is equal to the number rows of data.
##' 
##' @param train_params a GBMTrainParams object which specifies the
##' parameters used in growing decision trees.
##' 
##' @param var_monotone optional vector, the same length as the number
##' of predictors, indicating the relationship each variable has with
##' the outcome.  It have a monotone increasing (+1) or decreasing
##' (-1) or an arbitrary relationship.
##' 
##' @param var_names a vector of strings of containing the names of
##' the predictor variables.
##' 
##' @param cv_folds a positive integer specifying the number of folds
##' to be used in cross-validation of the gbm fit. If cv_folds > 1
##' then cross-validation is performed; the default of cv_folds is 1.
##' 
##' @param cv_class_stratify a bool specifying whether or not to
##' stratify via response outcome. Currently only applies to
##' "Bernoulli" distribution and defaults to false.
##' 
##' @param fold_id An optional vector of values identifying what fold
##' each observation is in. If supplied, cv_folds can be
##' missing. Note: Multiple rows of the same observation must have the
##' same fold_id.
##' 
##' @param keep_gbm_data a bool specifying whether or not the gbm_data
##' object created in this method should be stored in the results.
##' 
##' @param par_details Details of the parallelization to use in the
##' core algorithm (\code{\link{gbmParallel}}).
##'     
##' @param is_verbose if TRUE, gbmt will print out progress and
##' performance of the fit.
##' 
##' @return a \code{GBMFit} object.
##'
##' @examples
##' ## create some data
##' N <- 1000
##' X1 <- runif(N)
##' X2 <- runif(N)
##' X3 <- factor(sample(letters[1:4],N,replace=TRUE))
##' mu <- c(-1,0,1,2)[as.numeric(X3)]
##' 
##' p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
##' Y <- rbinom(N,1,p)
##' 
##' # random weights if you want to experiment with them
##' w <- rexp(N)
##' w <- N*w/sum(w)
##' 
##' data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
##' 
##' \dontrun{
##' train_params <-
##'      training_params(num_trees = 3000,
##'                      shrinkage = 0.001,
##'                      bag_fraction = 0.5,
##'                      num_train = N/2,
##'                      id=seq_len(nrow(data)),
##'                      min_num_obs_in_node = 10,
##'                      interaction_depth = 3,
##'                      num_features = 3)
##' }
##'
##' train_params <-
##'      training_params(num_trees = 100,
##'                      shrinkage = 0.001,
##'                      bag_fraction = 0.5,
##'                      num_train = N/2,
##'                      id=seq_len(nrow(data)),
##'                      min_num_obs_in_node = 10,
##'                      interaction_depth = 3,
##'                      num_features = 3)
##'  
##' # fit initial model
##' gbm1 <- gbmt(Y~X1+X2+X3,                # formula
##'             data=data,                 # dataset
##'             weights=w,
##'             var_monotone=c(0,0,0),     # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
##'             distribution=gbm_dist("Bernoulli"),
##'             train_params = train_params,
##'             cv_folds=5,                # do 5-fold cross-validation
##'             is_verbose = FALSE)           # don't print progress
##' 
##' # plot the performance
##' best.iter.oob <- gbmt_performance(gbm1,method="OOB")  # returns out-of-bag estimated best number of trees
##' plot(best.iter.oob)
##' print(best.iter.oob)
##' best.iter.cv <- gbmt_performance(gbm1,method="cv")   # returns 5-fold cv estimate of best number of trees
##' plot(best.iter.cv)
##' print(best.iter.cv)
##' best.iter.test <- gbmt_performance(gbm1,method="test") # returns test set estimate of best number of trees
##' plot(best.iter.cv)
##' print(best.iter.test)
##' 
##' best.iter <- best.iter.test
##' 
##' # plot variable influence
##' summary(gbm1,num_trees=1)         # based on the first tree
##' summary(gbm1,num_trees=best.iter) # based on the estimated best number of trees
##' 
##' # create marginal plots
##' # plot variable X1,X2,X3 after "best" iterations
##' par(mfrow=c(1,3))
##' plot(gbm1,1,best.iter)
##' plot(gbm1,2,best.iter)
##' plot(gbm1,3,best.iter)
##' par(mfrow=c(1,1))
##' plot(gbm1,1:2,best.iter) # contour plot of variables 1 and 2 after "best" number iterations
##' plot(gbm1,2:3,best.iter) # lattice plot of variables 2 and 3 after "best" number iterations
##' 
##' # 3-way plot
##' plot(gbm1,1:3,best.iter)
##' 
##' # print the first and last trees
##' print(pretty_gbm_tree(gbm1,1))
##' print(pretty_gbm_tree(gbm1, gbm1$params$num_trees))
##' @export
gbmt <- function(formula,
                 distribution=gbm_dist("Gaussian"),
                 data,
                 weights=rep(1, nrow(data)),
                 offset=rep(0, nrow(data)),
                 train_params=training_params(num_trees=2000,
                     interaction_depth=3,
                     min_num_obs_in_node=10,
                     shrinkage=0.001,
                     bag_fraction=0.5,
                     id=seq_len(nrow(data)),
                     num_train=round(0.5 * nrow(data)),
                     num_features=ncol(data)-1),
                 var_monotone=NULL,
                 var_names=NULL,
                 cv_folds=1,
                 cv_class_stratify=FALSE,
                 fold_id=NULL,
                 keep_gbm_data=FALSE,
                 par_details=getOption('gbm.parallel'),
                 is_verbose=FALSE) {
  # Extract the model
  the_call <- match.call()
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
  
  # If missing guess distribution
  if(missing(distribution)) {
    distribution <- gbm_dist(guess_distribution(y))
  }
  # Update distribution according to groups
  distribution <- determine_groups(data, y, distribution)
  
  # Update weights according to groupings
  weights <- weight_group_consistency(distribution, weights)
  
  # Update number of training rows based off of groups
  train_params <- update_num_train_groups(train_params, distribution)
  
  # Call gbmt.fit 
  gbm_fit_obj <- gbmt_fit(x, y, distribution, weights, offset,
                         train_params, as.character(formula[[2]]),
                         var_monotone, var_names, keep_gbm_data, cv_folds,
                         cv_class_stratify, fold_id, par_details, is_verbose)

  # Wrap up extra pieces 
  gbm_fit_obj$model <- m
  gbm_fit_obj$Terms <- Terms
  gbm_fit_obj$call <- the_call
  gbm_fit_obj$is_verbose <- is_verbose
 
  return(gbm_fit_obj)
}
