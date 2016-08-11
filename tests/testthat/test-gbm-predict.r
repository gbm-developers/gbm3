####################
# Author: James Hickey
#
# Series of test to check the behaviour of the predict S3 method.
#
####################

context("Test input checking - predict.GBMFit")
test_that("Error thrown if type is not 'link' or 'response'", {
  # Given a fit and data
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
  
  # Make prediction
  set.seed(2)
  # make some new data
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=TRUE))
  X4 <- ordered(sample(letters[1:6],N,replace=TRUE))
  X5 <- factor(sample(letters[1:3],N,replace=TRUE))
  X6 <- 3*runif(N)
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  # Actual underlying signal - Check how close to this we are
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  data2 <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # When calling predict with type not 'link' or 'response'
  # Then an error will be thrown
  expect_error(predict(fit, data2, length(fit$trees), type='wrong'), "type must be either 'link' or 'response'")
})
test_that("Error thrown if data is missing", {
  # Given a fit and data
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
  
  # Make prediction
  set.seed(2)
  # make some new data
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=TRUE))
  X4 <- ordered(sample(letters[1:6],N,replace=TRUE))
  X5 <- factor(sample(letters[1:3],N,replace=TRUE))
  X6 <- 3*runif(N)
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  # Actual underlying signal - Check how close to this we are
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  data2 <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # When calling predict without data
  # Then an error will be thrown
  expect_error(predict(fit, num_trees=length(fit$trees)), "new_data must be provided as a data frame")
})
test_that("Error thrown if data is not in data.frame", {
  # Given a fit and data
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
  
  # Make prediction
  set.seed(2)
  # make some new data
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=TRUE))
  X4 <- ordered(sample(letters[1:6],N,replace=TRUE))
  X5 <- factor(sample(letters[1:3],N,replace=TRUE))
  X6 <- 3*runif(N)
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  # Actual underlying signal - Check how close to this we are
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  data2 <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # When calling predict with data not in data.frame
  # Then an error will be thrown
  expect_error(predict(fit, new_data=as.list(data2), num_trees=length(fit$trees)), 
               "new_data must be provided as a data frame")
})
test_that("Error thrown if num_trees not provided", {
  # Given a fit and data
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
  
  # Make prediction
  set.seed(2)
  # make some new data
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=TRUE))
  X4 <- ordered(sample(letters[1:6],N,replace=TRUE))
  X5 <- factor(sample(letters[1:3],N,replace=TRUE))
  X6 <- 3*runif(N)
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  # Actual underlying signal - Check how close to this we are
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  data2 <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # When calling predict without number of trees specified
  # Then an error will be thrown
  expect_error(predict(fit, data2), "Number of trees to be used in prediction must be provided.")
})
test_that("Error thrown if num_trees is NULL or vector of 0 length", {
  # Given a fit and data
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
  
  # Make prediction
  set.seed(2)
  # make some new data
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=TRUE))
  X4 <- ordered(sample(letters[1:6],N,replace=TRUE))
  X5 <- factor(sample(letters[1:3],N,replace=TRUE))
  X6 <- 3*runif(N)
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  # Actual underlying signal - Check how close to this we are
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  data2 <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # When calling predict with number of trees NULL or length 0
  # Then an error will be thrown
  expect_error(predict(fit, data2, NULL), "num_trees cannot be NULL or a vector of zero length")
  expect_error(predict(fit, data2, c()), "num_trees cannot be NULL or a vector of zero length")
})
test_that("Error thrown if num_trees has element which is not an integer", {
  # Given a fit and data
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
  params <- training_params(num_trees=20000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE)
  
  # Make prediction
  set.seed(2)
  # make some new data
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=TRUE))
  X4 <- ordered(sample(letters[1:6],N,replace=TRUE))
  X5 <- factor(sample(letters[1:3],N,replace=TRUE))
  X6 <- 3*runif(N)
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  # Actual underlying signal - Check how close to this we are
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  data2 <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # When calling predict with number of trees containing an element which isn't an integer
  # Then an error will be thrown
  expect_error(predict(fit, data2, num_trees=c(1, -1.2)), "num_trees must be a vector of positive integers")
})
test_that("Warning thrown if offset was in original fit", {
  # Given a fit and data
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
  
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6 + offset(offset), data=data, distribution=dist, weights=w,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE)
  
  # Make prediction
  set.seed(2)
  # make some new data
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=TRUE))
  X4 <- ordered(sample(letters[1:6],N,replace=TRUE))
  X5 <- factor(sample(letters[1:3],N,replace=TRUE))
  X6 <- 3*runif(N)
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  # Actual underlying signal - Check how close to this we are
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  data2 <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # When calling predict with original fit containing the offset
  # Then a warning will be thrown
  expect_warning(predict(fit, data2, length(fit$trees)), 
                 "predict.GBMFit does not add the offset to the predicted values.")
})
test_that("Warning thrown if num_trees exceeds number in original fit", {
  # Given a fit and data
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
  
  # Make prediction
  set.seed(2)
  # make some new data
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=TRUE))
  X4 <- ordered(sample(letters[1:6],N,replace=TRUE))
  X5 <- factor(sample(letters[1:3],N,replace=TRUE))
  X6 <- 3*runif(N)
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  # Actual underlying signal - Check how close to this we are
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  data2 <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # When calling predict with number of trees exceeding number in original fit
  # Then a warning will be thrown
  expect_warning(predict(fit, data2, length(fit$trees)+1), 
                 "Number of trees exceeded number fit so far. Using ", paste(length(fit$trees),collapse=" "),".")
})

context("Test basic functionality of predict.GBMFit")
test_that("predict works if Terms defined in original fit", {
  # Given a fit and data
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
  
  # Make prediction
  set.seed(2)
  # make some new data
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=TRUE))
  X4 <- ordered(sample(letters[1:6],N,replace=TRUE))
  X5 <- factor(sample(letters[1:3],N,replace=TRUE))
  X6 <- 3*runif(N)
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  # Actual underlying signal - Check how close to this we are
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  data2 <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # When calling predict correctly
  # Then no error is thrown
  expect_error(predict(fit, data2, length(fit$trees)), NA)
})
test_that("type='response' scales predictions", {
  # Given a fit and data
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

  # make some new data
  set.seed(2)
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  data2 <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  
  # When calling predict with type = response
  preds <- predict(fit, data2, num_trees=length(fit$trees), type='response')
  preds_link <- predict(fit, data2, num_trees=length(fit$trees), type='link')
  
  # Then type=response predictions are scaled
  expect_equal(preds, adjust_pred_scale(preds_link, fit$distribution))
})
test_that("Output is matrix if length(num_trees) > 1", {
  # Given a fit and data
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
  
  # make some new data
  set.seed(2)
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  data2 <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  
  # When calling predict with length num_trees > 1
  preds <- predict(fit, data2, num_trees=c(50, 100))

  # Then output is a matrix
  expect_true(is.matrix(preds))
})
test_that("Output is vector if length(num_trees) == 1", {
  # Given a fit and data
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
  
  # make some new data
  set.seed(2)
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  data2 <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  
  # When calling predict with length num_trees = 1
  preds <- predict(fit, data2, num_trees=length(fit$trees), type='response')

  # Then output is a vector
  expect_true(is.vector(preds))
})
test_that("When num_trees specified exceeds total number in fit then number in fit used", {
  # Given a fit and data
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
  
  # make some new data
  set.seed(2)
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  data2 <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  
  # When calling predict with num_trees > length(fit$trees)
  # Then predictions are evaluated at num_trees = length(fit$trees)
  expect_equal(predict(fit, data2, num_trees=length(fit$trees) + 1),
               predict(fit, data2, num_trees=length(fit$trees)))
})