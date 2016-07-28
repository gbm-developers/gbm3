# Series of helper functions for plot.GBMFit

#### One variable helpers ####
#' @export get_ylabel_one_var
get_ylabel_one_var <- function(dist_obj) {
  UseMethod("get_ylabel_one_var", dist_obj)
}

get_ylabel_one_var.default <- function(dist_obj) {
  return("")
}

get_ylabel_one_var.BernoulliGBMDist <- function(dist_obj) {
  return("Predicted Probability")
}

get_ylabel_one_var.PairwiseGBMDist <- function(dist_obj) {
  return("Predicted Probability")
}

get_ylabel_one_var.PoissonGBMDist <- function(dist_obj) {
  return("Predicted Count")
}

#### Two variable helpers ####
select_two_var_plot <- function(f.factor, X, gbm_fit_obj, var_index, ...) {
  # Set plot identified
  which_plot <- sum(f.factor) + max(which(f.factor==TRUE))
  which_plot <- ifelse(length(which_plot)==0, 1, which_plot)
  which_plot <- toString(which_plot)
  
  # Call
  switch(which_plot,
         "1"=two_var_plot_no_factor(X, gbm_fit_obj, var_index, ...),
         "2"=two_var_plot_first_factor(X, gbm_fit_obj, var_index, ...),
         "3"=two_var_plot_second_factor(X, gbm_fit_obj, var_index, ...),
         "4"=two_var_plot_both_factor(X, gbm_fit_obj, var_index, ...))
}

two_var_plot_no_factor <- function(X, gbm_fit_obj, var_index, ...) {
  print(levelplot(y~X1*X2,data=X,
                  xlab=gbm_fit_obj$variables$var_names[var_index[1]],
                  ylab=gbm_fit_obj$variables$var.names[var_index[2]],...))
}

two_var_plot_first_factor <- function(X, gbm_fit_obj, var_index, ...) {
  print(xyplot(y~X2|X1,data=X,
               xlab=gbm_fit_obj$variables$var_names[var_index[2]],
               ylab=paste("f(", gbm_fit_obj$variables$var_names[var_index[1]],",", gbm_fit_obj$variables$var_names[var_index[2]],")",sep=""),
               type="l",
               panel = panel.xyplot,
               ...))
}

two_var_plot_second_factor <- function(X, gbm_fit_obj, var_index, ...) {
  print(xyplot(y~X1|X2,data=X,
               xlab=gbm_fit_obj$variables$var_names[var_index[1]],
               ylab=paste("f(",gbm_fit_obj$variables$var_names[var_index[1]],",",gbm_fit_obj$variables$var_names[var_index[2]],")",sep=""),
               type="l",
               panel = panel.xyplot,
               ...))
}

two_var_plot_both_factor <- function(X, gbm_fit_obj, var_index, ...) {
  print(stripplot(X1~y|X2,data=X,
                  xlab=gbm_fit_obj$variables$var_names[var_index[2]],
                  ylab=paste("f(",gbm_fit_obj$variables$var_names[var_index[1]],",",gbm_fit_obj$variables$var_names[var_index[2]],")",sep=""),
                  ...))
}

#### Three variable helpers ####
select_three_var_plot <- function(f.factor, X, gbm_fit_obj, var_index, ...) {
  which_plot <- toString(sum(f.factor))
  
  i <- order(f.factor)
  X.new <- X[,i]
  X.new$y <- X$y
  names(X.new) <- names(X)
  
  switch(which_plot,
         "0"=three_var_plot_no_factor(X.new, gbm_fit_obj, var_index, i,  ...),
         "1"=three_var_plot_one_factor(X.new, gbm_fit_obj, var_index, i, ...),
         "2"=three_var_plot_two_factor(X.new, gbm_fit_obj, var_index, i, ...),
         "3"=three_var_plot_three_factor(X.new, gbm_fit_obj, var_index, i, ...))
}

three_var_plot_no_factor <- function(X, gbm_fit_obj, var_index, select_index, ...) {

  X$X3 <- equal.count(X$X3)
  print(levelplot(y~X1*X2|X3,data=X,
                  xlab=gbm_fit_obj$variables$var_names[var_index[select_index[1]]],
                  ylab=gbm_fit_obj$variables$var_names[var_index[select_index[2]]],...))
  
}

three_var_plot_one_factor <- function(X, gbm_fit_obj, var_index, select_index, ...) {
  print(levelplot(y~X1*X2|X3,data=X,
                  xlab=gbm_fit_obj$variables$var_names[var_index[select_index[1]]],
                  ylab=gbm_fit_obj$variables$var_names[var_index[select_index[2]]],...))
}

three_var_plot_two_factor <- function(X, gbm_fit_obj, var_index, select_index, ...) {
  print(xyplot(y~X1|X2*X3,data=X,
               type="l",
               xlab=gbm_fit_obj$variables$var_names[var_index[select_index[1]]],
               ylab=paste("f(",paste(gbm_fit_obj$variables$var_names[var_index[1:3]],collapse=","),")",sep=""),
               panel = panel.xyplot,
               ...))
}

three_var_plot_three_factor <- function(X, gbm_fit_obj, var_index, select_index, ...) {
  print(stripplot(X1~y|X2*X3,data=X,
                  xlab=gbm_fit_obj$variables$var_names[var_index[select_index[1]]],
                  ylab=paste("f(",paste(gbm_fit_obj$variables$var_names[var_index[1:3]],collapse=","),")",sep=""),
                  ...))
}