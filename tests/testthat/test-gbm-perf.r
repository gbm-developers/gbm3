####################
# Author: James Hickey
#
# Series of tests to check gbm_perf
#
####################

context("gbm_perf input checking")
test_that("Error thrown if not passed GBMFit object", {
  # Given a fit object
  ## test Gaussian distribution gbm model
  set.seed(1)
  
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  X4 <- ordered(sample(letters[1:6],N,replace=T))
  X5 <- factor(sample(letters[1:3],N,replace=T))
  X6 <- 3*runif(N)
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  SNR <- 10 # signal-to-noise ratio
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  sigma <- sqrt(var(Y)/SNR)
  Y <- Y + rnorm(N,0,sigma)
  
  # create a bunch of missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  w <- rep(1,N)
  offset <- rep(0, N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  
  # Set up for new API
  params <- training_params(num_trees=200, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE)
  # When gbm_perf is called on fit object not of class GBMFit
  class(fit) <- "wrong"
  
  # Then an error is thrown
  expect_error(gbm_perf(fit, method="cv"))
})
test_that("Error thrown if plot_it is not a logical", {
  # Given a fit object
  ## test Gaussian distribution gbm model
  set.seed(1)
  
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  X4 <- ordered(sample(letters[1:6],N,replace=T))
  X5 <- factor(sample(letters[1:3],N,replace=T))
  X6 <- 3*runif(N)
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  SNR <- 10 # signal-to-noise ratio
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  sigma <- sqrt(var(Y)/SNR)
  Y <- Y + rnorm(N,0,sigma)
  
  # create a bunch of missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  w <- rep(1,N)
  offset <- rep(0, N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  
  # Set up for new API
  params <- training_params(num_trees=200, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE)
  
  # When calling gbm_perf with plot_it not a logical
  # Then an error is thrown
  expect_error(gbm_perf(fit, method="cv", plot_it=c(TRUE, FALSE)))
  expect_error(gbm_perf(fit, method="cv", plot_it=1.5))
  expect_error(gbm_perf(fit, method="cv", plot_it=""))
  expect_error(gbm_perf(fit, method="cv", plot_it=NaN))
  expect_error(gbm_perf(fit, method="cv", plot_it=Inf))
})
test_that("Error thrown if plot_it is NA", {
  # Given a fit object
  ## test Gaussian distribution gbm model
  set.seed(1)
  
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  X4 <- ordered(sample(letters[1:6],N,replace=T))
  X5 <- factor(sample(letters[1:3],N,replace=T))
  X6 <- 3*runif(N)
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  SNR <- 10 # signal-to-noise ratio
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  sigma <- sqrt(var(Y)/SNR)
  Y <- Y + rnorm(N,0,sigma)
  
  # create a bunch of missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  w <- rep(1,N)
  offset <- rep(0, N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  
  # Set up for new API
  params <- training_params(num_trees=200, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE)
  
  # When calling gbm_perf with plot_it=NA
  # Then error is thrown
  expect_error(gbm_perf(fit, method="cv", plot_it=NA))
})
test_that("Error thrown if method parameter missing", {
  # Given a fit object
  ## test Gaussian distribution gbm model
  set.seed(1)
  
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  X4 <- ordered(sample(letters[1:6],N,replace=T))
  X5 <- factor(sample(letters[1:3],N,replace=T))
  X6 <- 3*runif(N)
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  SNR <- 10 # signal-to-noise ratio
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  sigma <- sqrt(var(Y)/SNR)
  Y <- Y + rnorm(N,0,sigma)
  
  # create a bunch of missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  w <- rep(1,N)
  offset <- rep(0, N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  
  # Set up for new API
  params <- training_params(num_trees=200, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE)
  
  # When gbm_perf called with method missing
  # Then expect error
  expect_error(gbm_perf(fit))
})
test_that("Error thrown  if method not element of c('OOB', 'cv', 'test')", {
  # Given a fit object
  ## test Gaussian distribution gbm model
  set.seed(1)
  
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  X4 <- ordered(sample(letters[1:6],N,replace=T))
  X5 <- factor(sample(letters[1:3],N,replace=T))
  X6 <- 3*runif(N)
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  SNR <- 10 # signal-to-noise ratio
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  sigma <- sqrt(var(Y)/SNR)
  Y <- Y + rnorm(N,0,sigma)
  
  # create a bunch of missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  w <- rep(1,N)
  offset <- rep(0, N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  
  # Set up for new API
  params <- training_params(num_trees=200, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE)
  
  # When gbm_perf called with method that is not 'cv', 'test' or 'OOB'
  # Then an error is thrown
  expect_error(gbm_perf(fit, method="weird_metric"))
})
test_that("Message given if method is 'cv", {
  # Given a fit object
  ## test Gaussian distribution gbm model
  set.seed(1)
  
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  X4 <- ordered(sample(letters[1:6],N,replace=T))
  X5 <- factor(sample(letters[1:3],N,replace=T))
  X6 <- 3*runif(N)
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  SNR <- 10 # signal-to-noise ratio
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  sigma <- sqrt(var(Y)/SNR)
  Y <- Y + rnorm(N,0,sigma)
  
  # create a bunch of missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  w <- rep(1,N)
  offset <- rep(0, N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  
  # Set up for new API
  params <- training_params(num_trees=200, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE)
  
  # When gbm_perf is called with method 'cv'
  # Then a warning is thrown
  expect_message(gbm_perf(fit, method='cv'))
})
test_that("Message given if method is 'OOB'", {
  # Given a fit object
  ## test Gaussian distribution gbm model
  set.seed(1)
  
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  X4 <- ordered(sample(letters[1:6],N,replace=T))
  X5 <- factor(sample(letters[1:3],N,replace=T))
  X6 <- 3*runif(N)
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  SNR <- 10 # signal-to-noise ratio
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  sigma <- sqrt(var(Y)/SNR)
  Y <- Y + rnorm(N,0,sigma)
  
  # create a bunch of missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  w <- rep(1,N)
  offset <- rep(0, N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  
  # Set up for new API
  params <- training_params(num_trees=200, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE)
  
  # When gbm_perf is called with method 'OOB'
  # Then a warning is thrown
  expect_message(gbm_perf(fit, method="OOB"))
})

context("gbm_perf return")
test_that("gbm_perf returns correct best iteration for each method", {
  # Given a fit object and perf evaluated with each method
  ## test Gaussian distribution gbm model
  set.seed(1)
  
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  X4 <- ordered(sample(letters[1:6],N,replace=T))
  X5 <- factor(sample(letters[1:3],N,replace=T))
  X6 <- 3*runif(N)
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  SNR <- 10 # signal-to-noise ratio
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  sigma <- sqrt(var(Y)/SNR)
  Y <- Y + rnorm(N,0,sigma)
  
  # create a bunch of missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  w <- rep(1,N)
  offset <- rep(0, N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  
  # Set up for new API
  params <- training_params(num_trees=200, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE)
  
  best_iter_t <- which.min(fit$valid.error)
  best_iter_c <- which.min(fit$cv_error)
  x <- seq_len(fit$params$num_trees)
  smoother <- loess(fit$oobag.improve~x,
                    enp.target=min(max(4,length(x)/10),50))
  smoother$y <- smoother$fitted
  smoother$x <- x
  best_iter_oo <- smoother$x[which.min(-cumsum(smoother$y))]
  
  # When calling gbm_perf with 3 methods
  iter_t <- gbm_perf(fit, method="test")
  iter_c <- gbm_perf(fit, method="cv")
  iter_oo <- gbm_perf(fit, method="OOB")
  
  # Then correctly calculates best iterations
  expect_equal(iter_t, best_iter_t)
  expect_equal(iter_c, best_iter_c)
  expect_equal(iter_oo, best_iter_oo)
})