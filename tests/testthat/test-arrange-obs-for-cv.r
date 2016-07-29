####################
# Author: James Hickey
#
# Series of tests to check the extraction of observations
# into different folds
#
####################

context("Error checking")
test_that("split_and_join throws an error if gbm_data_obj is not GBMData object", {
  # Given gbm_data, a distribution obj, rows_in_fold and rows_in_training
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
  data <- data.frame(X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  
  # Set up for new API
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq_len(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  gdata <- gbm_data(data, Y, w, offset)
  cv_folds <- 5
  cv_groups <- create_cv_groups(gbm_data_obj = gdata, dist, params, cv_folds, FALSE, NULL)
  
  # Observations in the training set
  rows_in_training_set <- params$id %in% seq_len(params$num_train_rows)
  params$id <- params$id[rows_in_training_set]
  
  # Observations in cv_group - interested in first fold
  rows_in_fold <- params$id %in% seq_len(params$num_train_rows)[(cv_groups == 1)]
  
  # Extract relevent data - split into training and validation sets
  # Calculate new number of training rows
  params$num_train_rows <- length(which(cv_groups != 1))
  params$num_train <- length(unique(params$id[!rows_in_fold]))
  
  # When try to split_and_join not GBMData
  class(gdata) <- "wrong"
  
  # Then error is thrown
  expect_error(split_and_join(gdata, params, rows_in_training_set, rows_in_fold))
  
})
test_that("split_and_join throws an error if train_params is not GBMTrainParams object", {
  # Given gbm_data, a distribution obj, rows_in_fold and rows_in_training
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
  data <- data.frame(X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  
  # Set up for new API
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq_len(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  gdata <- gbm_data(data, Y, w, offset)
  cv_folds <- 5
  cv_groups <- create_cv_groups(gbm_data_obj = gdata, dist, params, cv_folds, FALSE, NULL)
  
  # Observations in the training set
  rows_in_training_set <- params$id %in% seq_len(params$num_train_rows)
  params$id <- params$id[rows_in_training_set]
  
  # Observations in cv_group - interested in first fold
  rows_in_fold <- params$id %in% seq_len(params$num_train_rows)[(cv_groups == 1)]
  
  # Extract relevent data - split into training and validation sets
  # Calculate new number of training rows
  params$num_train_rows <- length(which(cv_groups != 1))
  params$num_train <- length(unique(params$id[!rows_in_fold]))
  
  # When try to split_and_join not GBMTrainParams
  class(params) <- "wrong"
  
  # Then error is thrown
  expect_error(split_and_join(gdata, params, rows_in_training_set, rows_in_fold))
})
test_that("split_and_join throws an error if rows_in_training is not an atomic of logicals", {
  # Given gbm_data, a distribution obj, rows_in_fold and rows_in_training
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
  data <- data.frame(X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  
  # Set up for new API
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq_len(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  gdata <- gbm_data(data, Y, w, offset)
  cv_folds <- 5
  cv_groups <- create_cv_groups(gbm_data_obj = gdata, dist, params, cv_folds, FALSE, NULL)
  
  # Observations in the training set
  rows_in_training_set <- params$id %in% seq_len(params$num_train_rows)
  params$id <- params$id[rows_in_training_set]
  
  # Observations in cv_group - interested in first fold
  rows_in_fold <- params$id %in% seq_len(params$num_train_rows)[(cv_groups == 1)]
  
  # Extract relevent data - split into training and validation sets
  # Calculate new number of training rows
  params$num_train_rows <- length(which(cv_groups != 1))
  params$num_train <- length(unique(params$id[!rows_in_fold]))
  
  # When try to split_and_join with rows_in_training_set not an atomic of logicals
  # Then error is thrown
  expect_error(split_and_join(gdata, params, rows_in_training_set=c(1, 2), rows_in_fold))
  expect_error(split_and_join(gdata, params, rows_in_training_set="what", rows_in_fold))
  expect_error(split_and_join(gdata, params, rows_in_training_set= NaN, rows_in_fold))
  
})
test_that("split_and_join throws an error if rows_in_fold is not an atomic of logicals of length num_train", {
  # Given gbm_data, a distribution obj, rows_in_fold and rows_in_training
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
  data <- data.frame(X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  
  # Set up for new API
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq_len(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  gdata <- gbm_data(data, Y, w, offset)
  cv_folds <- 5
  cv_groups <- create_cv_groups(gbm_data_obj = gdata, dist, params, cv_folds, FALSE, NULL)
  
  # Observations in the training set
  rows_in_training_set <- params$id %in% seq_len(params$num_train_rows)
  params$id <- params$id[rows_in_training_set]
  
  # Observations in cv_group - interested in first fold
  rows_in_fold <- params$id %in% seq_len(params$num_train_rows)[(cv_groups == 1)]
  
  # Extract relevent data - split into training and validation sets
  # Calculate new number of training rows
  params$num_train_rows <- length(which(cv_groups != 1))
  params$num_train <- length(unique(params$id[!rows_in_fold]))
  
  # When try to split_and_join with rows_in_fold not an atomic of logicals of length num_train
  # Then error is thrown
  expect_error(split_and_join(gdata, params, rows_in_training_set, rows_in_fold =sample(c(TRUE, FALSE), params$num_train+1, replace=TRUE)))
  expect_error(split_and_join(gdata, params, rows_in_training_set, rows_in_fold=c(1, 2)))
  expect_error(split_and_join(gdata, params, rows_in_training_set, rows_in_fold=NaN))
  expect_error(split_and_join(gdata, params, rows_in_training_set, rows_in_fold="Where"))  
})
test_that("update_fold_dist_data throws an error if gbm_data_obj is not GBMData object", {
  # Given gbm_data, a distribution obj, rows_in_fold and rows_in_training
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
  data <- data.frame(X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  
  # Set up for new API
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq_len(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  gdata <- gbm_data(data, Y, w, offset)
  cv_folds <- 5
  cv_groups <- create_cv_groups(gbm_data_obj = gdata, dist, params, cv_folds, FALSE, NULL)
  
  # Observations in the training set
  rows_in_training_set <- params$id %in% seq_len(params$num_train_rows)
  params$id <- params$id[rows_in_training_set]
  
  # Observations in cv_group - interested in first fold
  rows_in_fold <- params$id %in% seq_len(params$num_train_rows)[(cv_groups == 1)]
  
  # Extract relevent data - split into training and validation sets
  # Calculate new number of training rows
  params$num_train_rows <- length(which(cv_groups != 1))
  params$num_train <- length(unique(params$id[!rows_in_fold]))
  
  # When try to update_fold_dist_data not GBMData
  class(gdata) <- "wrong"
  
  # Then error is thrown
  expect_error(update_fold_dist_data(dist, gdata, params, rows_in_training_set, rows_in_fold))
})
test_that("update_fold_dist_data throws an error if gbm_dist_obj is not GBMDist object", {
  # Given gbm_data, a distribution obj, rows_in_fold and rows_in_training
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
  data <- data.frame(X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  
  # Set up for new API
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq_len(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  gdata <- gbm_data(data, Y, w, offset)
  cv_folds <- 5
  cv_groups <- create_cv_groups(gbm_data_obj = gdata, dist, params, cv_folds, FALSE, NULL)
  
  # Observations in the training set
  rows_in_training_set <- params$id %in% seq_len(params$num_train_rows)
  params$id <- params$id[rows_in_training_set]
  
  # Observations in cv_group - interested in first fold
  rows_in_fold <- params$id %in% seq_len(params$num_train_rows)[(cv_groups == 1)]
  
  # Extract relevent data - split into training and validation sets
  # Calculate new number of training rows
  params$num_train_rows <- length(which(cv_groups != 1))
  params$num_train <- length(unique(params$id[!rows_in_fold]))
  
  # When try to update_fold_dist_data not GBMDist
  class(dist) <- "wrong"
  
  # Then error is thrown
  expect_error(update_fold_dist_data(dist, gdata, params, rows_in_training_set, rows_in_fold))
})
test_that("update_fold_dist_data throws an error if train_params is not GBMTrainParams object", {
  # Given gbm_data, a distribution obj, rows_in_fold and rows_in_training
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
  data <- data.frame(X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  
  # Set up for new API
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq_len(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  gdata <- gbm_data(data, Y, w, offset)
  cv_folds <- 5
  cv_groups <- create_cv_groups(gbm_data_obj = gdata, dist, params, cv_folds, FALSE, NULL)
  
  # Observations in the training set
  rows_in_training_set <- params$id %in% seq_len(params$num_train_rows)
  params$id <- params$id[rows_in_training_set]
  
  # Observations in cv_group - interested in first fold
  rows_in_fold <- params$id %in% seq_len(params$num_train_rows)[(cv_groups == 1)]
  
  # Extract relevent data - split into training and validation sets
  # Calculate new number of training rows
  params$num_train_rows <- length(which(cv_groups != 1))
  params$num_train <- length(unique(params$id[!rows_in_fold]))
  
  # When try to update_fold_dist_data not GBMTrainParams
  class(params) <- "wrong"
  
  # Then error is thrown
  expect_error(update_fold_dist_data(dist, gdata, params, rows_in_training_set, rows_in_fold))
})

context("Test splitting and joining functionality")
test_that("split_and_join returns a GBMData object", {
  # Given gbm_data, a distribution obj, rows_in_fold and rows_in_training
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
  data <- data.frame(X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  
  # Set up for new API
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq_len(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  gdata <- gbm_data(data, Y, w, offset)
  cv_folds <- 5
  cv_groups <- create_cv_groups(gbm_data_obj = gdata, dist, params, cv_folds, FALSE, NULL)
  
  # Observations in the training set
  rows_in_training_set <- params$id %in% seq_len(params$num_train_rows)
  params$id <- params$id[rows_in_training_set]
  
  # Observations in cv_group - interested in first fold
  rows_in_fold <- params$id %in% seq_len(params$num_train_rows)[(cv_groups == 1)]

  # Extract relevent data - split into training and validation sets
  # Calculate new number of training rows
  params$num_train_rows <- length(which(cv_groups != 1))
  params$num_train <- length(unique(params$id[!rows_in_fold]))
  
  # When split_and_join called
  returned_obj <- split_and_join(gdata, params, rows_in_training_set, rows_in_fold)
  
  # Then returned object is a GBMData object
  expect_error(check_if_gbm_data(returned_obj), NA)
})
test_that("split_and_join modifies gbm_data_obj so data in validation fold is at the 'end' of each field", {
  # Given gbm_data, a distribution obj, rows_in_fold and rows_in_training
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
  data <- data.frame(X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  
  # Set up for new API
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq_len(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  gdata <- gbm_data(data, Y, w, offset)
  cv_folds <- 5
  cv_groups <- create_cv_groups(gbm_data_obj = gdata, dist, params, cv_folds, FALSE, NULL)
  
  # Observations in the training set
  rows_in_training_set <- params$id %in% seq_len(params$num_train_rows)
  params$id <- params$id[rows_in_training_set]
  
  # Observations in cv_group - interested in first fold
  rows_in_fold <- params$id %in% seq_len(params$num_train_rows)[(cv_groups == 1)]
  
  # Extract relevent data - split into training and validation sets
  # Calculate new number of training rows
  params$num_train_rows <- length(which(cv_groups != 1))
  params$num_train <- length(unique(params$id[!rows_in_fold]))
  
  # When split_and_join called
  returned_obj <- split_and_join(gdata, params, rows_in_training_set, rows_in_fold)
  
  # Then returned object is a GBMData object with validation fold data at end of rows
  expect_equal(returned_obj$x[rows_in_training_set, ,drop=FALSE][(N/2 - length(which(cv_groups==1)) + 1):(N/2), ,drop=FALSE], 
               data[rows_in_training_set, ,drop=FALSE][rows_in_fold, ,drop=FALSE])
  expect_equal(returned_obj$y[rows_in_training_set, ,drop=FALSE][(N/2 - length(which(cv_groups==1)) + 1):(N/2), ,drop=FALSE]$V1, 
              Y[rows_in_training_set][rows_in_fold])
  expect_equal(returned_obj$weights[rows_in_training_set][(N/2 - length(which(cv_groups==1)) + 1):(N/2)], 
               w[rows_in_training_set][rows_in_fold])
  expect_equal(returned_obj$offset[rows_in_training_set][(N/2 - length(which(cv_groups==1)) + 1):(N/2)], 
               offset[rows_in_training_set][rows_in_fold])
})
test_that("split_and_join updates x_order so as to only use data in training folds", {
  # Given gbm_data, a distribution obj, rows_in_fold and rows_in_training
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
  data <- data.frame(X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  
  # Set up for new API
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq_len(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  gdata <- gbm_data(data, Y, w, offset)
  cv_folds <- 5
  cv_groups <- create_cv_groups(gbm_data_obj = gdata, dist, params, cv_folds, FALSE, NULL)
  
  # Observations in the training set
  rows_in_training_set <- params$id %in% seq_len(params$num_train_rows)
  params$id <- params$id[rows_in_training_set]
  
  # Observations in cv_group - interested in first fold
  rows_in_fold <- params$id %in% seq_len(params$num_train_rows)[(cv_groups == 1)]
  
  # Extract relevent data - split into training and validation sets
  # Calculate new number of training rows
  params$num_train_rows <- length(which(cv_groups != 1))
  params$num_train <- length(unique(params$id[!rows_in_fold]))
  
  # When split_and_join called
  returned_obj <- split_and_join(gdata, params, rows_in_training_set, rows_in_fold)
  
  # Get gdata training folds data
  gdata_train <- gdata
  gdata_train$x <- as.data.frame(gdata_train$x[rows_in_training_set, ,drop=FALSE][!rows_in_fold, ,drop=FALSE])
  gdata_train$y <- as.data.frame(as.matrix(gdata_train$y)[rows_in_training_set, ,drop=FALSE][!rows_in_fold, ,drop=FALSE])
  gdata_train$offset <- gdata_train$offset[rows_in_training_set][!rows_in_fold]
  gdata_train$weights <- gdata_train$weights[rows_in_training_set][!rows_in_fold]
  x_order_with_training_folds <- predictor_order(gdata_train, params)$x_order
  
  # Then returned object is a GBMData object with x_order updated so as to only 
  expect_equal(as.matrix(returned_obj$x_order), as.matrix(x_order_with_training_folds))
})
test_that("update_fold_dist_data returns the original distribution object if NOT Pairwise or CoxPH", {
  # Given gbm_data, a distribution obj, rows_in_fold and rows_in_training
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
  data <- data.frame(X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  
  # Set up for new API
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq_len(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  gdata <- gbm_data(data, Y, w, offset)
  cv_folds <- 5
  cv_groups <- create_cv_groups(gbm_data_obj = gdata, dist, params, cv_folds, FALSE, NULL)
  
  # Observations in the training set
  rows_in_training_set <- params$id %in% seq_len(params$num_train_rows)
  params$id <- params$id[rows_in_training_set]
  
  # Observations in cv_group - interested in first fold
  rows_in_fold <- params$id %in% seq_len(params$num_train_rows)[(cv_groups == 1)]
  
  # Extract relevent data - split into training and validation sets
  # Calculate new number of training rows
  params$num_train_rows <- length(which(cv_groups != 1))
  params$num_train <- length(unique(params$id[!rows_in_fold]))
  
  # When try to update_fold_dist_data with distributions which are not CoxPH or Pairwise
  # Not assuming no stratification in cv_groups or fold_id so same for all dists
  # Then nothing changes for those distributions
  expect_equal(update_fold_dist_data(gbm_dist("AdaBoost"), gdata, params, rows_in_training_set, rows_in_fold), gbm_dist("AdaBoost"))
  expect_equal(update_fold_dist_data(gbm_dist("Bernoulli"), gdata, params, rows_in_training_set, rows_in_fold), gbm_dist("Bernoulli"))
  expect_equal(update_fold_dist_data(gbm_dist("Gamma"), gdata, params, rows_in_training_set, rows_in_fold), gbm_dist("Gamma"))
  expect_equal(update_fold_dist_data(gbm_dist("Gaussian"), gdata, params, rows_in_training_set, rows_in_fold), gbm_dist("Gaussian"))
  expect_equal(update_fold_dist_data(gbm_dist("Huberized"), gdata, params, rows_in_training_set, rows_in_fold), gbm_dist("Huberized"))
  expect_equal(update_fold_dist_data(gbm_dist("Laplace"), gdata, params, rows_in_training_set, rows_in_fold), gbm_dist("Laplace"))
  expect_equal(update_fold_dist_data(gbm_dist("Poisson"), gdata, params, rows_in_training_set, rows_in_fold), gbm_dist("Poisson"))
  expect_equal(update_fold_dist_data(gbm_dist("Quantile"), gdata, params, rows_in_training_set, rows_in_fold), gbm_dist("Quantile"))
  expect_equal(update_fold_dist_data(gbm_dist("TDist"), gdata, params, rows_in_training_set, rows_in_fold), gbm_dist("TDist"))
  expect_equal(update_fold_dist_data(gbm_dist("Tweedie"), gdata, params, rows_in_training_set, rows_in_fold), gbm_dist("Tweedie"))
})
test_that("update_fold_dist_data updates the dist objects strata and sorted fields if CoxPH", {
  # Given gbm_data, a distribution obj (CoxPH!!), rows_in_fold and rows_in_training
  # Require Surv to be available
  require(survival)
  
  # create some data
  set.seed(1)
  N <- 3000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  f <- 0.5*sin(3*X1 + 5*X2^2 + mu/10)
  tt.surv <- rexp(N,exp(f))
  tt.cens <- rexp(N,0.5)
  delta <- as.numeric(tt.surv <= tt.cens)
  tt <- apply(cbind(tt.surv,tt.cens),1,min)
  
  # throw in some missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  # random weights if you want to experiment with them
  w <- rep(1,N)
  data <- data.frame(tt=tt,delta=delta,X1=X1,X2=X2,X3=X3)
  
  # Put into new API
  params <- training_params(num_trees=3000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5, id=seq_len(nrow(data)), num_train=N/2, num_features=3)
  
  
  # Set up for new API
  dist <- gbm_dist("CoxPH")
  gdata <- gbm_data(data, Surv(tt, delta), w, offset=rep(0, N))
  cv_folds <- 5
  cv_groups <- create_cv_groups(gbm_data_obj = gdata, dist, params, cv_folds, FALSE, NULL)
  
  # Observations in the training set
  rows_in_training_set <- params$id %in% seq_len(params$num_train_rows)
  params$id <- params$id[rows_in_training_set]
  
  # Observations in cv_group - interested in first fold
  rows_in_fold <- params$id %in% seq_len(params$num_train_rows)[(cv_groups == 1)]
  
  # Extract relevent data - split into training and validation sets
  # Calculate new number of training rows
  params$num_train_rows <- length(which(cv_groups != 1))
  params$num_train <- length(unique(params$id[!rows_in_fold]))
  
  # When try to update_fold_dist_data with distributions 
  dist_update <- update_fold_dist_data(dist, gdata, params, rows_in_training_set, rows_in_fold)
  
  # Then strata and sorted are updated
  expect_true(length(dist_update$strata) != length(dist$strata))
  expect_false(all(is.na(dist_update$sorted)))
})
test_that("update_fold_dist_data updates the dist objects group if Pairwise", {
  # Given gbm_data, a distribution obj (Pairwise!! + groups), rows_in_fold and rows_in_training
  # create some data - DOESN'T NEED TO REFLECT DIST
  set.seed(1)
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
  X <- data.frame(X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  data <- gbm_data(X, Y, w, offset)
  
  # Put into new API
  params <- training_params(num_trees=3000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5, id=seq_len(nrow(X)), num_train=N/2, num_features=3)
  
  
  # Set up for new API
  dist <- gbm_dist("Pairwise")
  dist$group <- as.factor(sample(seq_len(5), N, replace=TRUE))
  cv_folds <- 5
  cv_groups <- create_cv_groups(gbm_data_obj = data, dist, params, cv_folds, FALSE, NULL)
 
  
  # DEFINE GROUPINGS
  # Observations in the training set
  rows_in_training_set <- params$id %in% seq_len(params$num_train_rows)
  params$id <- params$id[rows_in_training_set]
  
  # Observations in cv_group - interested in first fold
  rows_in_fold <- params$id %in% seq_len(params$num_train_rows)[(cv_groups == 1)]
  
  # Extract relevent data - split into training and validation sets
  # Calculate new number of training rows
  params$num_train_rows <- length(which(cv_groups != 1))
  params$num_train <- length(unique(params$id[!rows_in_fold]))
  
  # When try to update_fold_dist_data with distributions 
  dist_update <- update_fold_dist_data(dist, data, params, rows_in_training_set, rows_in_fold)
  
  # Then groupings are now updated
  expect_equal(dist_update$group,  c(dist$group[rows_in_training_set][!rows_in_fold],
                                     dist$group[rows_in_training_set][rows_in_fold]))
})
test_that("extract_obs_in_fold updates the training parameters so as to take into account the number of observations in the validation fold", {
  # Given gbm_data, a distribution obj (CoxPH!!), rows_in_fold and rows_in_training
  # Require Surv to be available
  require(survival)
  
  # create some data
  set.seed(1)
  N <- 3000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  f <- 0.5*sin(3*X1 + 5*X2^2 + mu/10)
  tt.surv <- rexp(N,exp(f))
  tt.cens <- rexp(N,0.5)
  delta <- as.numeric(tt.surv <= tt.cens)
  tt <- apply(cbind(tt.surv,tt.cens),1,min)
  
  # throw in some missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  # random weights if you want to experiment with them
  w <- rep(1,N)
  data <- data.frame(tt=tt,delta=delta,X1=X1,X2=X2,X3=X3)
  
  # Put into new API
  dist <- gbm_dist("CoxPH", prior_node_coeff_var=10)
  params <- training_params(num_trees=3000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=3)
  
  
  # Set up for new API
  dist <- gbm_dist("CoxPH")
  gdata <- gbm_data(data, Surv(tt, delta), w, offset=rep(0, N))
  cv_folds <- 5
  cv_groups <- create_cv_groups(gbm_data_obj = gdata, dist, params, cv_folds, FALSE, NULL)
  
  # When extract_obs_in_fold is called
  returned_obj <- extract_obs_in_fold(gdata, dist, params, cv_groups, fold_num=1) 
  
  # Then updates the num of training observations in train_params
  expect_equal(returned_obj$params$num_train, params$num_train * (cv_folds-1)/cv_folds)
})
test_that("extract_obs_in_fold returns a list of the updated gbm_data_obj, dist obj and training params", {
  # Given gbm_data, a distribution obj (CoxPH!!), rows_in_fold and rows_in_training
  # Require Surv to be available
  require(survival)
  
  # create some data
  set.seed(1)
  N <- 3000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  f <- 0.5*sin(3*X1 + 5*X2^2 + mu/10)
  tt.surv <- rexp(N,exp(f))
  tt.cens <- rexp(N,0.5)
  delta <- as.numeric(tt.surv <= tt.cens)
  tt <- apply(cbind(tt.surv,tt.cens),1,min)
  
  # throw in some missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  # random weights if you want to experiment with them
  w <- rep(1,N)
  data <- data.frame(X1=X1,X2=X2,X3=X3)
  
  # Put into new API
  params <- training_params(num_trees=3000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=3)
  
  
  # Set up for new API
  dist <- gbm_dist("CoxPH")
  gdata <- gbm_data(data, Surv(tt, delta), w, rep(0, N))
  cv_folds <- 5
  cv_groups <- create_cv_groups(gdata, dist, params, cv_folds, FALSE, NULL)
  
  # When extract_obs_in_fold is called
  returned_obj <- extract_obs_in_fold(gdata, dist, params, cv_groups, fold_num=1) 
    
  # Then returns a list with correct fields
  expect_true(is.list(returned_obj))
  expect_equal(names(returned_obj), c("data", "dist", "params"))
  expect_error(check_if_gbm_data(returned_obj$data), NA)
  expect_error(check_if_gbm_dist(returned_obj$dist), NA)
  expect_error(check_if_gbm_train_params(returned_obj$params), NA)
})



