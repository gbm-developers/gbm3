####################
# Author: James Hickey
#
# Series of tests to check print and summary methods
#
####################

context("Testing print.GBMFit helper functions")
test_that("print_perf_measures defaults to total number of iterations if train_fraction = 1 and fit is not cross-validated", {
  # Given a "correct" GBMFit object - train_fraction=1
  set.seed(1)
  
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  
  # random weights if you want to experiment with them
  w <- rexp(N)
  w <- N*w/sum(w)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  offset <- rep(0, N)
  
  # Set up for new API
  params <- training_params(num_trees=3000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5, id=seq(nrow(data)), num_train=N, num_features=3)
  dist <- gbm_dist("Bernoulli")
  fit <- gbmt(Y~X1+X2+X3, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0), keep_gbm_data=TRUE, cv_folds=1, is_verbose = FALSE)
  
  # When calculating print_perf_measures
  best_iter <- print_perf_measures(fit)
  
  # Then results equal to num_trees
  expect_equal(best_iter, params$num_trees)
  
})

test_that("print_perf_measures calculates the performance using test if train_fraction < 1", {
  # Given a "correct" GBMFit object - train_fraction < 1
  set.seed(1)
  
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  
  # random weights if you want to experiment with them
  w <- rexp(N)
  w <- N*w/sum(w)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  offset <- rep(0, N)
  
  # Set up for new API
  params <- training_params(num_trees=3000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=3)
  dist <- gbm_dist("Bernoulli")
  fit <- gbmt(Y~X1+X2+X3, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0), keep_gbm_data=TRUE, cv_folds=5, is_verbose = FALSE)
  
  # When calculating print_perf_measures
  best_iter <- print_perf_measures(fit)
  
  # Then results equal to gbm_perf with "test" method
  expect_equal(best_iter, gbm_perf(fit, method="test"))
})

test_that("print_perf_measures returns the performance using cv if fit is cross-validated and train_fraction=1", {
  # Given a "correct" GBMFit object - train_fraction=1 and cv_folds > 1
  set.seed(1)
  
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  
  # random weights if you want to experiment with them
  w <- rexp(N)
  w <- N*w/sum(w)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  offset <- rep(0, N)
  
  # Set up for new API
  params <- training_params(num_trees=3000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5, id=seq(nrow(data)), num_train=N, num_features=3)
  dist <- gbm_dist("Bernoulli")
  fit <- gbmt(Y~X1+X2+X3, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0), keep_gbm_data=TRUE, cv_folds=5, is_verbose = FALSE)
  
  # When calculating print_perf_measures
  best_iter <- print_perf_measures(fit)
  
  # Then results equal to gbm_perf with method="cv"
  expect_equal(best_iter, gbm_perf(fit, method="cv"))
})

test_that("print_iters_and_dist does not throw error when passed a GBMFit object", {
  # Given a "correct" GBMFit object
  set.seed(1)
  
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  
  # random weights if you want to experiment with them
  w <- rexp(N)
  w <- N*w/sum(w)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  offset <- rep(0, N)
  
  # Set up for new API
  params <- training_params(num_trees=3000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=3)
  dist <- gbm_dist("Bernoulli")
  fit <- gbmt(Y~X1+X2+X3, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0), keep_gbm_data=TRUE, cv_folds=5, is_verbose = FALSE)
  
  # When calculating print_iters_and_dist 
  # Then no error is thrown
  expect_error(print_iters_and_dist(fit), NA)
})

test_that("print_confusion_matrix does not throw error when passed a GBMFit object", {
  # Given a "correct" GBMFit object
  set.seed(1)
  
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  
  # random weights if you want to experiment with them
  w <- rexp(N)
  w <- N*w/sum(w)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  offset <- rep(0, N)
  
  # Set up for new API
  params <- training_params(num_trees=3000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=3)
  dist <- gbm_dist("Bernoulli")
  fit <- gbmt(Y~X1+X2+X3, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0), keep_gbm_data=TRUE, cv_folds=5, is_verbose = FALSE)
  
  # When calculating print_confusion_matrix
  # Then no error is thrown
  expect_error(print_confusion_matrix(fit), NA)
})

test_that("binary_response_conf_matrix does not throw error when called correctly", {
  # Given a "correct" GBMFit object - Bernoulli distribution
  set.seed(1)
  
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  
  # random weights if you want to experiment with them
  w <- rexp(N)
  w <- N*w/sum(w)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  offset <- rep(0, N)
  
  # Set up for new API
  params <- training_params(num_trees=3000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=3)
  dist <- gbm_dist("Bernoulli")
  fit <- gbmt(Y~X1+X2+X3, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0), keep_gbm_data=TRUE, cv_folds=5, is_verbose = FALSE)
  
  # When calculating binary_response_conf_matrix
  # Then no error is thrown
  expect_error(binary_response_conf_matrix(fit$gbm_data_obj$y, fit$cv_fitted), NA)
})

test_that("pseudo_r_squared does not throw error when passed correct inputs", {
  # Given a "correct" GBMFit object - using Gauss
  set.seed(1)
  
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  
  # random weights if you want to experiment with them
  w <- rexp(N)
  w <- N*w/sum(w)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  offset <- rep(0, N)
  
  # Set up for new API
  params <- training_params(num_trees=3000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=3)
  dist <- gbm_dist("Bernoulli")
  fit <- gbmt(Y~X1+X2+X3, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0), keep_gbm_data=TRUE, cv_folds=5, is_verbose = FALSE)
  
  # When calculating pseudo_r_squared 
  # Then no error is thrown
  expect_error(pseudo_r_squared(fit$gbm_data_obj$y, fit$cv_fitted, fit$distribution$name), NA)
})

test_that("pseudo_r_squared is same for all dists except Gaussian", {
  # Given a "correct" GBMFit object - NOT Gaussian
  set.seed(1)
  
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  
  # random weights if you want to experiment with them
  w <- rexp(N)
  w <- N*w/sum(w)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  offset <- rep(0, N)
  
  # Set up for new API
  params <- training_params(num_trees=3000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=3)
  dist <- gbm_dist("Bernoulli")
  fit <- gbmt(Y~X1+X2+X3, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0), keep_gbm_data=TRUE, cv_folds=5, is_verbose = FALSE)
  
  # When calculating pseudo_r_squared for Bernoulli
  r_sq <- pseudo_r_squared(fit$gbm_data_obj$y, fit$cv_fitted, fit$distribution$name)
  
  # Then this value is same for all dists except Gaussian
  expect_equal(pseudo_r_squared(fit$gbm_data_obj$y, fit$cv_fitted, "AdaBoost"), r_sq)
  expect_equal(pseudo_r_squared(fit$gbm_data_obj$y, fit$cv_fitted, "Bernoulli"), r_sq)
  expect_equal(pseudo_r_squared(fit$gbm_data_obj$y, fit$cv_fitted, "CoxPH"), r_sq)
  expect_equal(pseudo_r_squared(fit$gbm_data_obj$y, fit$cv_fitted, "Gamma"), r_sq)
  expect_true(pseudo_r_squared(fit$gbm_data_obj$y, fit$cv_fitted, "Gaussian") != r_sq)
  expect_equal(pseudo_r_squared(fit$gbm_data_obj$y, fit$cv_fitted, "Huberized"), r_sq)
  expect_equal(pseudo_r_squared(fit$gbm_data_obj$y, fit$cv_fitted, "Laplace"), r_sq)
  expect_equal(pseudo_r_squared(fit$gbm_data_obj$y, fit$cv_fitted, "Pairwise"), r_sq)
  expect_equal(pseudo_r_squared(fit$gbm_data_obj$y, fit$cv_fitted, "Poisson"), r_sq)
  expect_equal(pseudo_r_squared(fit$gbm_data_obj$y, fit$cv_fitted, "Quantile"), r_sq)
  expect_equal(pseudo_r_squared(fit$gbm_data_obj$y, fit$cv_fitted, "TDist"), r_sq)
  expect_equal(pseudo_r_squared(fit$gbm_data_obj$y, fit$cv_fitted, "Tweedie"), r_sq)
})


context("Test Input Error checking - printing")
test_that("print_perf_measures throws error if passed an object other than GBMFit class", {
  # Given a "correct" GBMFit objects
  set.seed(1)
  
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  
  # random weights if you want to experiment with them
  w <- rexp(N)
  w <- N*w/sum(w)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  offset <- rep(0, N)
  
  # Set up for new API
  params <- training_params(num_trees=3000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=3)
  dist <- gbm_dist("Bernoulli")
  fit <- gbmt(Y~X1+X2+X3, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0), keep_gbm_data=TRUE, cv_folds=5, is_verbose = FALSE)
  
  
  # When  print_perf_measures is called on fit which is missing class GBMFit
  class(fit) <- ""
  
  # Then an error is thrown
  expect_error(print_perf_measures(fit))
})

test_that("print_iters_and_dist throws error if passed an object other than GBMFit class", {
  # Given a "correct" GBMFit objects
  set.seed(1)
  
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  
  # random weights if you want to experiment with them
  w <- rexp(N)
  w <- N*w/sum(w)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  offset <- rep(0, N)
  
  # Set up for new API
  params <- training_params(num_trees=3000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=3)
  dist <- gbm_dist("Bernoulli")
  fit <- gbmt(Y~X1+X2+X3, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0), keep_gbm_data=TRUE, cv_folds=5, is_verbose = FALSE)
  
  
  # When  print_iters_and_dist is called on fit which is missing class GBMFit
  class(fit) <- ""
  
  # Then an error is thrown
  expect_error(print_iters_and_dist(fit))
})

test_that("print_confusion_matrix throws error if passed an object other than GBMFit class", {
  # Given a "correct" GBMFit objects
  set.seed(1)
  
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  
  # random weights if you want to experiment with them
  w <- rexp(N)
  w <- N*w/sum(w)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  offset <- rep(0, N)
  
  # Set up for new API
  params <- training_params(num_trees=3000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=3)
  dist <- gbm_dist("Bernoulli")
  fit <- gbmt(Y~X1+X2+X3, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0), keep_gbm_data=TRUE, cv_folds=5, is_verbose = FALSE)
  
  
  # When  print_confustion_matrix is called on fit which is missing class GBMFit
  class(fit) <- ""
  
  # Then an error is thrown
  expect_error(print_confusion_matrix(fit))
})

test_that("summary.GBMFit throws error if cBars is not a whole number > 1", {
  # Given a "correct" GBMFit objects
  set.seed(1)
  
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  
  # random weights if you want to experiment with them
  w <- rexp(N)
  w <- N*w/sum(w)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  offset <- rep(0, N)
  
  # Set up for new API
  params <- training_params(num_trees=3000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=3)
  dist <- gbm_dist("Bernoulli")
  fit <- gbmt(Y~X1+X2+X3, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0), keep_gbm_data=TRUE, cv_folds=5, is_verbose = FALSE)
  
  
  # When summary method is called with cBars is not a whole number > 1
  # Then an error is thrown
  expect_error(summary(fit, cBars=c(TRUE, FALSE)))
  expect_error(summary(fit, cBars=1.3))
  expect_error(summary(fit, cBars="What"))
  expect_error(summary(fit, cBars=Inf))
  expect_error(summary(fit, cBars=NULL))
})

test_that("summary.GBMFit throws error if num_trees is not a whole number > 1", {
  # Given a "correct" GBMFit objects
  set.seed(1)
  
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  
  # random weights if you want to experiment with them
  w <- rexp(N)
  w <- N*w/sum(w)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  offset <- rep(0, N)
  
  # Set up for new API
  params <- training_params(num_trees=3000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=3)
  dist <- gbm_dist("Bernoulli")
  fit <- gbmt(Y~X1+X2+X3, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0), keep_gbm_data=TRUE, cv_folds=5, is_verbose = FALSE)
  
  
  # When summary method is called with num_trees is not a whole number > 1
  # Then an error is thrown
  expect_error(summary(fit, num_trees=c(TRUE, FALSE)))
  expect_error(summary(fit, num_trees=1.3))
  expect_error(summary(fit, num_trees="What"))
  expect_error(summary(fit, num_trees=Inf))
  expect_error(summary(fit, num_trees=NULL))
})

test_that("summary.GBMFit throws error if plot_it is not a logical", {
  # Given a "correct" GBMFit objects
  set.seed(1)
  
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  
  # random weights if you want to experiment with them
  w <- rexp(N)
  w <- N*w/sum(w)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  offset <- rep(0, N)
  
  # Set up for new API
  params <- training_params(num_trees=3000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=3)
  dist <- gbm_dist("Bernoulli")
  fit <- gbmt(Y~X1+X2+X3, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0), keep_gbm_data=TRUE, cv_folds=5, is_verbose = FALSE)
  
  
  # When summary method is called with plot_it not a logical
  # Then an error is thrown
  expect_error(summary(fit, plot_it=c(TRUE, FALSE)))
  expect_error(summary(fit, plot_it=1))
  expect_error(summary(fit, plot_it="What"))
  expect_error(summary(fit, plot_it=Inf))
  expect_error(summary(fit, plot_it=NULL))
})

test_that("summary.GBMFit throws error if normalize is not a logical", {
  # Given a "correct" GBMFit objects
  set.seed(1)
  
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  
  # random weights if you want to experiment with them
  w <- rexp(N)
  w <- N*w/sum(w)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  offset <- rep(0, N)
  
  # Set up for new API
  params <- training_params(num_trees=3000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=3)
  dist <- gbm_dist("Bernoulli")
  fit <- gbmt(Y~X1+X2+X3, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0), keep_gbm_data=TRUE, cv_folds=5, is_verbose = FALSE)
 
  
  # When summary method is called with normalize not a logical
  # Then an error is thrown
  expect_error(summary(fit, normalize=c(TRUE, FALSE)))
  expect_error(summary(fit, normalize=1))
  expect_error(summary(fit, normalize="What"))
  expect_error(summary(fit, normalize=Inf))
  expect_error(summary(fit, normalize=NULL))
})

context("Test print.GBMFit")
test_that("Print method on GBMFit object runs without error", {
  # Given some "correct" GBMFit objects
  set.seed(1)
  
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  
  # random weights if you want to experiment with them
  w <- rexp(N)
  w <- N*w/sum(w)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  offset <- rep(0, N)
  
  # Set up for new API
  params <- training_params(num_trees=3000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=3)
  dist <- gbm_dist("Bernoulli")
  fit <- gbmt(Y~X1+X2+X3, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0), keep_gbm_data=TRUE, cv_folds=5, is_verbose = FALSE)
  fit2 <- gbmt(Y~X1+X2+X3, data=data, distribution=dist, weights=w, offset=offset,
               train_params=params, var_monotone=c(0, 0, 0), keep_gbm_data=TRUE, cv_folds=1, is_verbose = FALSE)
  
  # When print method is called
  # Then runs without error
  expect_error(print(fit), NA)
  expect_error(print(fit2), NA)
})
context("Test summary.GBMFit")
test_that("Summary method on GBMFit object runs without error", {
  # Given some "correct" GBMFit objects
  set.seed(1)
  
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  
  # random weights if you want to experiment with them
  w <- rexp(N)
  w <- N*w/sum(w)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  offset <- rep(0, N)
  
  # Set up for new API
  params <- training_params(num_trees=3000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=3)
  dist <- gbm_dist("Bernoulli")
  fit <- gbmt(Y~X1+X2+X3, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0), keep_gbm_data=TRUE, cv_folds=5, is_verbose = FALSE)
  fit2 <- gbmt(Y~X1+X2+X3, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0), keep_gbm_data=TRUE, cv_folds=1, is_verbose = FALSE)

  # When summary method is called
  # Then runs without error
  expect_error(summary(fit), NA)
  expect_error(summary(fit2), NA)
})

test_that("Summary method returns data.frame of variables and relative influences", {
  # Given some "correct" GBMFit objects
  set.seed(1)
  
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  
  # random weights if you want to experiment with them
  w <- rexp(N)
  w <- N*w/sum(w)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  offset <- rep(0, N)
  
  # Set up for new API
  params <- training_params(num_trees=3000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=3)
  dist <- gbm_dist("Bernoulli")
  fit <- gbmt(Y~X1+X2+X3, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0), keep_gbm_data=TRUE, cv_folds=5, is_verbose = FALSE)
  fit2 <- gbmt(Y~X1+X2+X3, data=data, distribution=dist, weights=w, offset=offset,
               train_params=params, var_monotone=c(0, 0, 0), keep_gbm_data=TRUE, cv_folds=1, is_verbose = FALSE)
  
  # When summary method is called
  summary_fit_1 <- summary(fit)
  summary_fit_2 <- summary(fit2)
  
  # Then summary method returns data.frame of variables and relative influences
  # with variables ordered 
  expect_true(is.data.frame(summary_fit_1))
  expect_true(is.data.frame(summary_fit_2))
  expect_equal(names(summary_fit_1), c("var", "rel_inf"))
  expect_equal(names(summary_fit_2), c("var", "rel_inf"))
})

test_that("Summary method returns variables and relative influence ordered by relative influence in descending order", {
  # Given a "correct" GBMFit objects
  set.seed(1)
  
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  
  # random weights if you want to experiment with them
  w <- rexp(N)
  w <- N*w/sum(w)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  offset <- rep(0, N)
  
  # Set up for new API
  params <- training_params(num_trees=3000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=3)
  dist <- gbm_dist("Bernoulli")
  fit <- gbmt(Y~X1+X2+X3, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0), keep_gbm_data=TRUE, cv_folds=5, is_verbose = FALSE)
  
  # When summary method is called
  summary_fit_1 <- summary(fit)
  
  # Then the variables are ordered in terms of descending relative influence
  rel_inf <- relative_influence(fit, gbm_perf(fit, method="cv"))
  expect_equal(as.factor(fit$variables$var_names[order(-rel_inf)]), summary_fit_1$var)
})