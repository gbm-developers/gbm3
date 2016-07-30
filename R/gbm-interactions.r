#' Estimate the strength of interaction effects
#' 
#' Computes Friedman's H-statistic to assess the strength of variable
#' interactions.
#' 
#' \code{interact.GBMFit} computes Friedman's H-statistic to assess the relative
#' strength of interaction effects in non-linear models. H is on the scale of
#' [0-1] with higher values indicating larger interaction effects. To connect
#' to a more familiar measure, if \eqn{x_1} and \eqn{x_2} are uncorrelated
#' covariates with mean 0 and variance 1 and the model is of the form
#' \deqn{y=\beta_0+\beta_1x_1+\beta_2x_2+\beta_3x_3} then
#' \deqn{H=\frac{\beta_3}{\sqrt{\beta_1^2+\beta_2^2+\beta_3^2}}}
#' 
#' Note that if the main effects are weak, the estimated H will be unstable.
#' For example, if (in the case of a two-way interaction) neither main effect
#' is in the selected model (relative influence is zero), the result will be
#' 0/0. Also, with weak main effects, rounding errors can result in values of H
#' > 1 which are not possible.
#' 
#' @usage interact(gbm_fit_obj, data, var_indices=1, num_trees=gbm_fit_obj$params$num_trees)
#' 
#' @param gbm_fit_obj a \code{GBMFit} object fitted using a call to \code{\link{gbmt}}.
#' 
#' @param data the dataset used to construct \code{gbm_fit_obj}. If the original dataset
#' is large, a random subsample may be used to accelerate the computation.
#' 
#' @param var_indices a vector of indices or the names of the variables for compute
#' the interaction effect. If using indices, the variables are indexed in the
#' same order that they appear in the initial \code{gbmt} formula.
#' 
#' @param num_trees the number of trees used to generate the plot. Only the first
#' \code{num_trees} trees will be used.
#' 
#' @return Returns the value of \code{H}.
#' 
#' @author Greg Ridgeway \email{gregridgeway@@gmail.com}
#' @seealso \code{\link{gbmt}}
#' @references J.H. Friedman and B.E. Popescu (2005). \dQuote{Predictive
#' Learning via Rule Ensembles.} Section 8.1
#' @keywords methods
#' @export 
#'

interact <- function(gbm_fit_obj, data, var_indices=1, num_trees=gbm_fit_obj$params$num_trees) {
  UseMethod("interact", gbm_fit_obj)
}

#' @name interact
#' @export
interact.GBMFit <- function(gbm_fit_obj, data, var_indices=1, num_trees = gbm_fit_obj$params$num_trees){
  # Initial input checks and set up
  if(!is.data.frame(data) && !is.matrix(data)) {
    stop("data argument should be a data.frame or matrix")  
  }
  if(!is.atomic(var_indices) ||
     !(all(var_indices == as.integer(var_indices)) || all(var_indices == as.character(var_indices)))) {
    stop("Observation ids must be a vector of integers or characters")
  }
  if (gbm_fit_obj$params$interaction_depth < length(var_indices)){
    stop("interaction_depth < length(variables_indices): too low in model call")
  }
  var_indices <- check_and_set_variables_indices(gbm_fit_obj, var_indices)
  num_trees <- check_and_set_num_trees(gbm_fit_obj, num_trees)

  # Convert factors to appropriate numerics
  for(var in var_indices) {
    if(is.factor(data[, gbm_fit_obj$variables$var_names[var]]))
      data[, gbm_fit_obj$variables$var_names[var]] <-
        as.numeric(data[, gbm_fit_obj$variables$var_names[var]])-1
  }
  
  # Generate a list with all combinations of variables
  all_combinations_vars <- apply(expand.grid(rep(list(c(FALSE,TRUE)), length(var_indices)))[-1,], 1,
             function(x) as.numeric(which(x)))
  
  # Compute predictions and "parity" for all variable combinations  
  preds_for_comb_vars <- compute_preds_for_all_var_combinations(gbm_fit_obj, all_combinations_vars, var_indices, num_trees)
  
  # Compute H-statistic
  # Set to prediction with all variables
  H_stat_squared <- preds_for_comb_vars[[length(all_combinations_vars)]]$preds
  
  # Loop over other combinations and see what variables have been excluded
  # Add to predictions for all variables with correct sign
  for(vars in seq_len((length(all_combinations_vars)-1))){
    i1 <- apply(preds_for_comb_vars[[length(all_combinations_vars)]]$data[, all_combinations_vars[[vars]], drop=FALSE],
                1, paste, collapse="\r")
    i2 <- apply(preds_for_comb_vars[[vars]]$data, 1, paste,collapse="\r")
    i <- match(i1, i2)
    
    H_stat_squared <- H_stat_squared + with(preds_for_comb_vars[[vars]], sign*preds[i,])
  }
  
  # The H-statistic squared is given by sum over variables predictions with all included 
  # minus the partial dependence on a variable + partial dependence excluding a variable
  # This sum is normalized by the sum of the prediction with no variables excluded
  weights <- matrix(preds_for_comb_vars[[length(all_combinations_vars)]]$num_levels_factors, ncol=1)
  sum_preds_no_exclusion <- matrix(preds_for_comb_vars[[length(all_combinations_vars)]]$preds^2, ncol=1, byrow=FALSE)
  
  numerator <- apply(H_stat_squared^2, 2, weighted.mean, w = weights, na.rm = TRUE)
  denominator <- apply(sum_preds_no_exclusion, 2, weighted.mean, w = weights, na.rm = TRUE)
  H_stat_squared <- numerator / denominator
  
  # If H > 1, rounding and tiny main effects have messed things up
  H_stat_squared[H_stat_squared > 1] <- NaN
  
  return(sqrt(H_stat_squared))
}


#### Helper Functions - Not to be used outside of this interact function ####
check_and_set_num_trees <- function(gbm_fit_obj, num_trees) {
  if(length(num_trees) > 1) {
    warning("length num_trees > 1: using first element")
    num_trees <- num_trees[1]
  }
  check_if_natural_number(num_trees)
  if (num_trees > gbm_fit_obj$params$num_trees) {
    warning(paste("num_trees exceeds the number of trees in the model, ",
                  gbm_fit_obj$params$num_trees,". Using ", gbm_fit_obj$params$num_trees, " trees.", sep = ""))
    num_trees <- gbm_fit_obj$params$num_trees
  }
  
  return(num_trees)
}

check_and_set_variables_indices <- function(gbm_fit_obj, variables_indices) {
  # Match up variable_indices to var_names - convert characters
  if (all(is.character(variables_indices))){
    i <- match(variables_indices, gbm_fit_obj$variables$var_names)
    if (any(is.na(i))) {
      stop("Variables given are not used in gbm model fit: ", variables_indices[is.na(i)])
    }
    else {
      variables_indices <- i
    }
  }
  
  if ((min(variables_indices) < 1) || (max(variables_indices) > length(gbm_fit_obj$variables$var_names))) {
    warning("variables_indices must be between 1 and ", length(gbm_fit_obj$variables$var_names))
  }
  
  return(variables_indices)
}

table_of_unique_values <- function(data, variables_indices) {
  unique_vars <- unique(data[, variables_indices,drop=FALSE])
  unique_vars$num_levels_factors <- table(factor(apply(data[, variables_indices,drop=FALSE],1,paste,collapse="\r"),
                      levels=apply(unique_vars, 1,paste,collapse="\r")))
  return(unique_vars)
}

compute_preds_for_all_var_combinations <- function(gbm_fit_obj, all_combinations_vars, variables_indices, num_trees) {
  preds_for_comb_vars <- vector("list", length(all_combinations_vars))
  for(vars in seq_along(all_combinations_vars)) {
    # Get data for combination
    preds_for_comb_vars[[vars]]$data <- data.frame(table_of_unique_values(data, 
                                                                          gbm_fit_obj$variables$var_names[variables_indices[all_combinations_vars[[vars]]]]))
    preds_for_comb_vars[[vars]]$num_levels_factors <- as.numeric(preds_for_comb_vars[[vars]]$data$num_levels_factors)
    preds_for_comb_vars[[vars]]$data$num_levels_factors <- NULL
    
    # Make predictions using the current combination of variables
    preds_for_comb_vars[[vars]]$preds <- .Call("gbm_plot",
                                               X = data.matrix(preds_for_comb_vars[[vars]]$data),
                                               i.var = as.integer(variables_indices[all_combinations_vars[[vars]]] - 1),
                                               n.trees = as.integer(num_trees),
                                               initF = as.double(gbm_fit_obj$initF),
                                               trees = gbm_fit_obj$trees,
                                               c.splits = gbm_fit_obj$c.splits,
                                               var.type = as.integer(gbm_fit_obj$variables$var_type),
                                               PACKAGE = "gbm")
    
    # Convert predictions to flat matrix
    preds_for_comb_vars[[vars]]$preds <- matrix(preds_for_comb_vars[[vars]]$preds, ncol=1, byrow=FALSE)
    
    # Centre the predictions
    preds_for_comb_vars[[vars]]$preds <- apply(preds_for_comb_vars[[vars]]$preds, 2, function(x, w){
      x - weighted.mean(x, w, na.rm=TRUE)
    }, w=preds_for_comb_vars[[vars]]$num_levels_factors)
    
    # precompute the sign of these terms to appear in H - statistic
    # if same "parity" return 1, else -1
    preds_for_comb_vars[[vars]]$sign <- ifelse((length(all_combinations_vars[[vars]]) %% 2) == 
                                                        (length(variables_indices) %% 2), 1, -1)
  }
  
  return(preds_for_comb_vars)
}
