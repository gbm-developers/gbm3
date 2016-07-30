####################
# Author: James Hickey
#
# Series of tests for gbm_more
#
####################

context("Checking inputs to gbm_more")
test_that("Error thrown if gbm_fit_obj is not of class GBMFit", {
  # Given a GBMFit from data
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
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE)
  
  # When stripped of GBMFit class and gbm_more is called
  class(fit) <- "wrong"
  
  # Then an error is thrown
  expect_error(gbm_more(fit))
})
test_that("Error thrown if num_new_trees is not a natural number > 1", {
  # Given an appropriate GBMFit object
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
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE)
  # When gbm_more is called with num_new_trees not a natural number > 1
  # Then an error is thrown
  expect_error(gbm_more(fit, num_new_trees = c(1, 200)))
  expect_error(gbm_more(fit, num_new_trees = FALSE))
  expect_error(gbm_more(fit, num_new_trees = NA))
  expect_error(gbm_more(fit, num_new_trees = Inf))
  expect_error(gbm_more(fit, num_new_trees = NaN))
  expect_error(gbm_more(fit, num_new_trees = "Eh?"))
  expect_error(gbm_more(fit, num_new_trees = 10.3))
})
test_that("Error thrown if is_verbose is not a logical", {
  # Given an appropriate GBMFit object
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
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE)
  # When gbm_more is called with is_verbose not a logical
  # Then an error is thrown
  expect_error(gbm_more(fit, is_verbose = c(1, 200)))
  expect_error(gbm_more(fit, is_verbose = 100))
  expect_error(gbm_more(fit, is_verbose =  NA))
  expect_error(gbm_more(fit, is_verbose =  Inf))
  expect_error(gbm_more(fit, is_verbose =  NaN))
  expect_error(gbm_more(fit, is_verbose = "Eh?"))
  expect_error(gbm_more(fit, is_verbose = 10.3))
})
test_that("Error thrown if data not provided or gbm_data_obj not stored in gbm_fit_obj", {
  # Given an appropriate GBMFit object - DATA NOT KEPT
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
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=keep_data, cv_folds=10, is_verbose=FALSE)
  # When gbm_more is called without providing new data
  # Then an error is thrown
  expect_error(gbm_more(fit))
})
test_that("Error thrown if data has different number of rows than the original number of rows in response", {
  # Given an appropriate GBMFit object (DATA NOT KEPT) 
  # and new data with more rows than before
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
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=keep_data, cv_folds=10, is_verbose=FALSE)
  
  # make some new data - more than last time
  N <- N + 1
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=TRUE))
  X4 <- ordered(sample(letters[1:6],N,replace=TRUE))
  X5 <- factor(sample(letters[1:3],N,replace=TRUE))
  X6 <- 3*runif(N)
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  # Actual underlying signal - Check how close to this we are
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  new_data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # When gbm_more is called with this new_data
  # Then an error is thrown
  expect_error(gbm_more(fit, data = new_data))
})

context("Check gbm_more properties")
test_that("Output of gbm_more is GBMFit object", {
  # Given an appropriate GBMFit object (DATA KEPT) 
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
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=keep_data, cv_folds=10, is_verbose=FALSE)
  
 
  
  # When gbm_more is called 
  fit_2 <- gbm_more(fit)
  
  # Then returned object is a GBMFit object
  expect_true(class(fit_2) %in% "GBMFit")
})
test_that("gbm_more does not update cv properties of GBMFit output", {
  # Given an appropriate GBMFit object (DATA KEPT) 
  # produced from cv
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
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=keep_data, cv_folds=10, is_verbose=FALSE)
  
  
  
  # When gbm_more is called 
  fit_2 <- gbm_more(fit)
  
  # Then returned objects cv fields have not changed
  expect_equal(fit_2$cv_error, fit$cv_error)
  expect_equal(fit_2$cv_fold, fit$cv_folds)
  expect_equal(fit_2$cv_fitted, fit$cv_fitted)
})
test_that("gbm_more fits additional trees - reflected in length of fit fields", {
  # Given an appropriate GBMFit object (DATA KEPT) 
  # with number of new trees
  keep_data <- TRUE
  num_trees <- 150
  
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
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=keep_data, cv_folds=10, is_verbose=FALSE)
  
  
  
  # When gbm_more is called 
  fit_2 <- gbm_more(fit, num_new_trees = num_trees)
  
  # Then errors and trees have these additional fitted trees to the object
  expect_equal(length(fit_2$train.error) - length(fit$train.error), num_trees)
  expect_equal(length(fit_2$valid.error) - length(fit$valid.error), num_trees)
  expect_equal(length(fit_2$oobag.improve) - length(fit$oobag.improve), num_trees)
  expect_equal(length(fit_2$trees) - length(fit$trees), num_trees)
})
test_that("Can run gbm_more with new data", {
  # Given an appropriate GBMFit object (DATA NOT KEPT) 
  # and new data same number of rows as before
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
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=keep_data, cv_folds=10, is_verbose=FALSE)
  
  # make some new data - more than last time
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=TRUE))
  X4 <- ordered(sample(letters[1:6],N,replace=TRUE))
  X5 <- factor(sample(letters[1:3],N,replace=TRUE))
  X6 <- 3*runif(N)
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  # Actual underlying signal - Check how close to this we are
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  new_data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # When calling gbm_more with new data
  # Then runs successfully
  expect_error(gbm_more(fit, data=new_data), NA)
})
test_that("Output of gbm_more DOES NOT contain constructed gbm_data_obj if run with new data", {
  # Given an appropriate GBMFit object (DATA NOT KEPT) 
  # and new data same number of rows as before
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
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=keep_data, cv_folds=10, is_verbose=FALSE)
  
  # make some new data - more than last time
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=TRUE))
  X4 <- ordered(sample(letters[1:6],N,replace=TRUE))
  X5 <- factor(sample(letters[1:3],N,replace=TRUE))
  X6 <- 3*runif(N)
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  # Actual underlying signal - Check how close to this we are
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  new_data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # When calling gbm_more with new data
  fit_more <- gbm_more(fit, data=new_data)
  
  # Then output dpes not store new_data object
  expect_true(is.null(fit_more$gbm_data_obj))
})
test_that("If run with new data and weights are length 0 then they are set to 1 for each data row", {
  # Given an appropriate GBMFit object (DATA KEPT) 
  # and new data same number of rows as before
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
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=keep_data, cv_folds=10, is_verbose=FALSE)
  
  # make some new data - more than last time
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=TRUE))
  X4 <- ordered(sample(letters[1:6],N,replace=TRUE))
  X5 <- factor(sample(letters[1:3],N,replace=TRUE))
  X6 <- 3*runif(N)
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  # Actual underlying signal - Check how close to this we are
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  new_data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # When calling gbm_more with new data and weights length 0
  fit_more <- gbm_more(fit, data=new_data, weights = NULL)
  
  # Then output weights are vector of 1s
  expect_equal(fit_more$gbm_data_obj$weights, rep(1, N))
})
test_that("If run with new data and offset is length 0 then they are set to 0 for each data row", {
  # Given an appropriate GBMFit object (DATA KEPT) 
  # and new data same number of rows as before
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
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=keep_data, cv_folds=10, is_verbose=FALSE)
  
  # make some new data - more than last time
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=TRUE))
  X4 <- ordered(sample(letters[1:6],N,replace=TRUE))
  X5 <- factor(sample(letters[1:3],N,replace=TRUE))
  X6 <- 3*runif(N)
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  # Actual underlying signal - Check how close to this we are
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  new_data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # When calling gbm_more with new data and offset length 0
  fit_more <- gbm_more(fit, data=new_data, offset = NULL)
  
  # Then output offset are vector of 0s
  expect_equal(fit_more$gbm_data_obj$offset, rep(0, N))
})