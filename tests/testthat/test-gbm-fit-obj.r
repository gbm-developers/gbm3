####################
# Author: James Hickey
#
# Series of test to validate the GBMFit object
#
####################

context("Test GBMFit definition")
test_that("Calling gbmt creates a GBMFit object", {
  # Given Data from examples
  ## Based on example in R package
  
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
  params <- training_params(num_trees=20, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  
  # When gbmt called
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE)
  
  
  # Then fit has correct class
  expect_equal(class(fit), "GBMFit")
})
test_that("GBMFit object has correct fields", {
  # Given Data from examples
  ## Based on example in R package
  
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
  params <- training_params(num_trees=20, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  
  # When gbmt called
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE)
  
  
  # Then fit has names
  expect_equal(names(fit), c("initF", "fit", "train.error", "valid.error", "oobag.improve", "trees", "c.splits", "distribution", "params"       
               ,"variables", "cv_folds", "cv_error", "cv_fitted", "gbm_data_obj",
               "par_details", "response_name", "model", "Terms", "call", "is_verbose"))
})
test_that("GBMFit contains correct GBM objects", {
  # Given Data from examples
  ## Based on example in R package
  
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
  params <- training_params(num_trees=20, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  
  # When gbmt called
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE)
  
  
  # Then fit has correct GBM objects
  expect_error(check_if_gbm_data(fit$gbm_data_obj), NA)
  expect_error(check_if_gbm_dist(fit$distribution), NA)
  expect_error(check_if_gbm_train_params(fit$params), NA)
  expect_error(check_if_gbm_var_container(fit$variables), NA)  
})
test_that("GBMFit does not contain gbm_data_obj if keep_gbm_data=FALSE", {
  # Given Data from examples
  ## Based on example in R package
  keep_data <- FALSE
  
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
  params <- training_params(num_trees=20, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  
  # When gbmt called
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=keep_data, cv_folds=10, is_verbose=FALSE)
  
  
  # Then fit does not store gbm data from gbmt
  expect_true(is.null(fit$gbm_data_obj))
})
test_that("call to gbm and gbmt produces object with same class and fields", {
  # Given Data from examples
  ## Based on example in R package
  keep_data <- TRUE
  
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
  params <- training_params(num_trees=20, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  
  # When gbmt called and gbm called
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=keep_data, cv_folds=10, is_verbose=FALSE)
  fit_2 <- gbm(Y~X1+X2+X3+X4+X5+X6, data=data, weights=w, offset=offset, n.trees = 2000, cv.folds = 10, keep.data = keep_data, verbose = FALSE)

  
  # Then same object and fields produced
  expect_equal(class(fit), class(fit_2))
  expect_equal(names(fit), names(fit_2))
})