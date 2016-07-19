#' Print model summary
#' 
#' Display basic information about a \code{GBMFit} object.
#' 
#' Prints some information about the model object. In particular, this method
#' prints the call to \code{gbm2}, the type of loss function that was used,
#' and the total number of iterations.
#' 
#' If cross-validation was performed, the 'best' number of trees as estimated
#' by cross-validation error is displayed. If a test set was used, the 'best'
#' number of trees as estimated by the test set error is displayed.
#' 
#' The number of available predictors, and the number of those having non-zero
#' influence on predictions is given (which might be interesting in data mining
#' applications).
#' 
#' If bernoulli or adaboost was used, the confusion matrix and
#' prediction accuracy are printed (objects being allocated to the class with
#' highest probability for bernoulli). These classifications
#' are performed using the cross-validation fitted values.
#' 
#' If the 'distribution' was specified as gaussian, laplace, quantile or
#' t-distribution, a summary of the residuals is displayed.  The residuals are
#' the cross-validation residuals. Also, a pseudo R-squared value is displayed.
#' For Gaussian response, this is 1 - sum(r*r) / sum(z*z) where z = y -
#' mean(y). For the other distributions, this is 1 - (median(abs(r)) /
#' mad(y))^2, following the suggestion of Rousseeuw and Leroy (equation 3.11).
#' Note that this definition of a robust R-squared is contentious.
#' 
#' @usage print(gbm_fit_obj, ...)
#' 
#' @param gbm_fit_obj is a fitted generalized boosting object of class \code{GBMFit}.
#' 
#' @param \dots arguments passed to \code{print.default}.
#' 
#' @author Harry Southworth, Daniel Edwards
#' @seealso \code{\link{gbm2}}
#' @references P. J. Rousseeuw and A. M. Leroy, Robust Regression and Outlier
#' Detection, Wiley, 1987 (2003).
#' @keywords models nonlinear survival nonparametric
#' @export print.GBMFit
#' 

print.GBMFit <- function(gbm_fit_obj, ... ){
  #  Print out number of iterations and distribution used
  print_iters_and_dist(gbm_fit_obj)
  
  # Print out performance measures
  best_iter <- print_perf_measures(gbm_fit_obj)
  
  # Print out relative influence of variables
  ri <- relative_influence(gbm_fit_obj, num_trees=best_it)
  cat( "There were", length(gbm_fit_obj$variables$var_names), "predictors of which",
       sum(ri > 0), "had non-zero influence.\n" )
  
  # CV confusion matrix and pseudo-R-squared
  print_confusion_mat(gbm_fit_obj)
  
  return(invisible())
}


#### Helper Functions ####
print_iters_and_dist <- function(gbm_fit_obj) {
  # If Pairwise dist extract metric and max_rank
  # else distribution details is just the name
  check_if_gbm_fit(gbm_fit_obj)
  if(gbm_fit_obj$distribution$name == "Pairwise") {
    if (!is.null(gbm_fit_obj$distribution$max.rank) && (gbm_fit_obj$distribution$max.rank > 0)) {
      distribution_details <- sprintf("pairwise (metric=%s, max.rank=%d)", gbm_fit_obj$distribution$metric,
                                      gbm_fit_obj$distribution$max.rank)
    } else {
      distribution_details <- sprintf("pairwise (metric=%s)", gbm_fit_obj$distribution$metric)
    }
  } else {
    distribution_details <- gbm_fit_obj$distribution$name
  }
  cat( paste( "A gradient boosted model with", distribution_details, "loss function.\n" ))
  cat( paste( length(gbm_fit_obj$train.error), "iterations were performed.\n" ) )
}

print_perf_measures <- function(gbm_fit_obj) {
  # Calculate the best number of iterations - returns test set if 
  # possible
  check_if_gbm_fit(gbm_fit_obj)
  
  # Set default answer - final iteration
  best_iter <- length(gbm_fit_obj$train_error)
  
  # CV best iteration 
  if (!is.null(gbm_fit_obj$cv_error)) {
    best_iter <- gbm_perf(gbm_fit_obj, plot_it = FALSE, method="cv" )
    cat(paste("The best cross-validation iteration was ", best_iter, ".\n", sep = "" ))
  }
  
  # Test set best iteration
  if (gbm_fit_obj$params$train_fraction < 1 ) {
    best_iter <- gbm_perf(gbm_fit_obj, plot_it = FALSE, method="test" )
    cat( paste("The best test-set iteration was ", best_iter, ".\n", sep = "" ) )
  }
  
  return(best_iter)
}

print_confusion_matrix <- function(gbm_fit_obj) {
  # This prints the confusion matrix based on the
  # cross validated fit - if no cv fit is present - break
  if(is.null(gbm_fit_obj$cv_fitted)) return(invisible())
  
  # If data was not kept than can't calculate confusion matrix
  check_if_gbm_fit(gbm_fit_obj)
  if(is.null(gbm_fit_obj$gbm_data_obj)) return(invisible())
  
  # Print off confusion matrix or pseudo R^2
  if (gbm_fit_obj$distribution$name %in% c("Bernoulli", "AdaBoost", "Huberized")) {
    binary_response_conf_matrix(gbm_fit_obj$gbm_data_obj$y, gbm_fit_obj$cv_fitted)
  } else if (gbm_fit_obj$distribution$name %in% c("Gaussian", "Laplace", "Quantile", "TDist")) {
    pseudo_r_squared(gbm_fit_obj$gbm_data_obj$y, gbm_fit_obj$cv_fitted, gbm_fit_obj$distribution$name)
  }
}

binary_response_conf_matrix <- function(response, cv_fit) {
  # Function to print out the confusion matrix for a binary response 
  # cross-validated generalized boosted model
  p <- 1 / (1 + exp(-cv_fit))
  p <- ifelse(p < .5, 0, 1)
  
  conf_mat <- matrix(table(c(response + 2 * p , 0:3)), ncol=2)
  conf_mat <- conf.mat - 1
  pred_acc <- round(100 * sum(diag(conf_mat)) / sum(conf_mat),2)
  
  conf_mat <- cbind(conf_mat,  round(100*diag(conf_mat)/rowSums(conf_mat),2))
  dimnames(conf_mat) <- list(c("0","1"), c("0", "1", "Pred. Acc."))
  
  cat("\nCross-validation confusion matrix:\n")
  print(conf_mat)
  cat("\nCross-validation prediction Accuracy = ", pred_acc, "%\n", sep = "")
}

pseudo_r_squared <- function(response, cv_fit, dist_name) {
  # Calculate residuals
  cv_resids <- response - gbm_fit_obj$cv_fitted
  
  cat("\nSummary of cross-validation residuals:\n" )
  print(quantile(cv_resids))
  cat("\n")
  
  # Do pseudo R^2
  if (dist_name == "Gaussian"){
    yadj <- response - mean(response[, 1])
    R2 <- 1 - sum(r^2)/sum(yadj^2)
    cat("Cross-validation pseudo R-squared: ", signif(R2, 3), "\n")
  }
  else { # Rousseeuw & Leroy, page 44
    R2 <- 1 - (median(abs(cv_resids)) / mad(response[, 1]))^2
    cat("Cross-validation robust pseudo R-squared: ", signif(R2, 3), "\n")
  }
}

