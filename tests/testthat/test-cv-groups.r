####################
# Author: James Hickey
#
# Series of tests to check the creation of cv groups
#
####################

context("Testing input checking - create_cv_groups")
test_that("create_cv_groups throws an error if not passed a GBMData object", {
  # Given inputs
  # Dist_Obj
  dist <- gbm_dist()
  
  # create some data
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
  
  # Params Obj
  params <- training_params(num_train = N)
  
  # Other CV Parameters
  cv_class_stratify <- FALSE
  cv_folds <- 2
  fold_id <- NULL
  
  # When data is stripped of GBMData class
  class(data) <- "wrong"
  
  # Then error thrown when trying to calculate CV groups
  expect_error(create_cv_groups(data, dist, params, cv_folds, cv_class_stratify, fold_id))
})
test_that("create_cv_groups throws an error if not passed a GBMDist object", {
  # Given inputs
  # Dist_Obj
  dist <- gbm_dist()
  
  # create some data
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
  
  # Params Obj
  params <- training_params(num_train = N)
  
  # Other CV Parameters
  cv_class_stratify <- FALSE
  cv_folds <- 2
  fold_id <- NULL
  
  # When dist is stripped of GBMDist class
  class(dist) <- ""
  
  # Then error thrown when trying to calculate CV groups
  expect_error(create_cv_groups(data, dist, params, cv_folds, cv_class_stratify, fold_id))
})
test_that("create_cv_groups throws an error if not passed a GBMTrainParams object", {
  # Given inputs
  # Dist_Obj
  dist <- gbm_dist()
  
  # create some data
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
  
  # Params Obj
  params <- training_params(num_train = N)
  
  # Other CV Parameters
  cv_class_stratify <- FALSE
  cv_folds <- 2
  fold_id <- NULL
  
  # When params is stripped of GBMTrainParams class
  class(params) <- ""
  
  # Then error thrown when trying to calculate CV groups
  expect_error(create_cv_groups(data, dist, params, cv_folds, cv_class_stratify, fold_id))
})
test_that("create_cv_groups throws an error if not cv_folds is not a whole number >= 1", {
  # Given inputs - cv_folds not a whole number >= 1
  # Dist_Obj
  dist <- gbm_dist()
  
  # create some data
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
  
  # Params Obj
  params <- training_params(num_train = N)
  
  # Other CV Parameters
  cv_class_stratify <- FALSE
  fold_id <- NULL
  cv_folds_1 <- 2.1
  cv_folds_2 <- -1
  cv_folds_3 <- "What?"
  cv_folds_4 <- Inf
  cv_folds_5 <- NaN
  cv_folds_6 <- 0
  
  # When creating cv groups
  # Then error thrown when trying to calculate CV groups
  expect_error(create_cv_groups(data, dist, params, cv_folds_1, cv_class_stratify, fold_id))
  expect_error(create_cv_groups(data, dist, params, cv_folds_2, cv_class_stratify, fold_id))
  expect_error(create_cv_groups(data, dist, params, cv_folds_3, cv_class_stratify, fold_id))
  expect_error(create_cv_groups(data, dist, params, cv_folds_4, cv_class_stratify, fold_id))
  expect_error(create_cv_groups(data, dist, params, cv_folds_5, cv_class_stratify, fold_id))
  expect_error(create_cv_groups(data, dist, params, cv_folds_6, cv_class_stratify, fold_id))
})
test_that("create_cv_groups throws an error if not cv_class_stratify is not a logical", {
  # Given inputs - cv_class_stratify is not a logical!
  # Dist_Obj
  dist <- gbm_dist()
  
  # create some data
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
  
  # Params Obj
  params <- training_params(num_train = N)
  
  # Other CV Parameters
  fold_id <- NULL
  cv_folds <- 2
  cv_class_stratify_1 <- -1
  cv_class_stratify_2 <- "What?"
  cv_class_stratify_3 <- Inf
  cv_class_stratify_4 <- NaN
  cv_class_stratify_5 <- 0
  cv_class_stratify_6 <- c(TRUE, FALSE)
  
  # When creating cv groups
  # Then error thrown when trying to calculate CV groups
  expect_error(create_cv_groups(data, dist, params, cv_folds, cv_class_stratify_1, fold_id))
  expect_error(create_cv_groups(data, dist, params, cv_folds, cv_class_stratify_2, fold_id))
  expect_error(create_cv_groups(data, dist, params, cv_folds, cv_class_stratify_3, fold_id))
  expect_error(create_cv_groups(data, dist, params, cv_folds, cv_class_stratify_4, fold_id))
  expect_error(create_cv_groups(data, dist, params, cv_folds, cv_class_stratify_5, fold_id))
  expect_error(create_cv_groups(data, dist, params, cv_folds, cv_class_stratify_6, fold_id))
})

context("Testing cv groupings")
test_that("cv groups is an atomic vector of integers whose length is equal to the number of training observations", {
  # Given inputs
  # Dist_Obj
  dist <- gbm_dist()
  
  # create some data
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
  
  # Params Obj
  params <- training_params(num_train = N)
  
  # Other CV Parameters
  cv_class_stratify <- FALSE
  cv_folds <- 2
  fold_id <- NULL
  
  # When creating cv_groups
  cv_groups <- create_cv_groups(data, dist, params, cv_folds, cv_class_stratify, fold_id)
  
  # Then cv_groups is an atomic vector of integers with length equal to N
  expect_true(any(cv_groups == as.integer(cv_groups)))
  expect_true(is.atomic(cv_groups))
  expect_true(length(cv_groups)== params$num_train)  
})
test_that("max integer in cv groups produced is = cv_folds parameter", {
  # Given inputs
  # Dist_Obj
  dist <- gbm_dist()
  
  # create some data
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
  
  # Params Obj
  params <- training_params(num_train = N)
  
  # Other CV Parameters
  cv_class_stratify <- FALSE
  cv_folds <- 2
  fold_id <- NULL
  
  # When creating cv_groups
  cv_groups <- create_cv_groups(data, dist, params, cv_folds, cv_class_stratify, fold_id)
  
  # Then cv_groups max fold number == cv_folds
  expect_true(max(cv_groups) == cv_folds)
})
test_that("cv groups is equal to fold_id if not NULL - distributions other than Bernoulli or Pairwise", {
  # Given inputs - default distributions other than Bernoulli
  # Dist_Objs
  dist_1 <- gbm_dist("AdaBoost")
  dist_2 <- gbm_dist("CoxPH")
  dist_3 <- gbm_dist("Gamma")
  dist_4 <- gbm_dist("Gaussian")
  dist_5 <- gbm_dist("Huberized")
  dist_6 <- gbm_dist("Laplace")
  dist_7 <- gbm_dist("Poisson")
  dist_8 <- gbm_dist("Quantile")
  dist_9 <- gbm_dist("TDist")
  dist_10 <- gbm_dist("Tweedie")
  
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
  
  # Params Obj
  params <- training_params(num_train = N)
  
  # Other CV Parameters - NB: fold_id not NULL
  cv_class_stratify <- FALSE
  cv_folds <- 2
  fold_id <- seq_len(params$num_train)
  
  # When creating cv_groups for distributions other than Bernoulli
  cv_groups_1 <- create_cv_groups(data, dist_1, params, cv_folds, cv_class_stratify, fold_id)
  cv_groups_2 <- create_cv_groups(data, dist_2, params, cv_folds, cv_class_stratify, fold_id)
  cv_groups_3 <- create_cv_groups(data, dist_3, params, cv_folds, cv_class_stratify, fold_id)
  cv_groups_4 <- create_cv_groups(data, dist_4, params, cv_folds, cv_class_stratify, fold_id)
  cv_groups_5 <- create_cv_groups(data, dist_5, params, cv_folds, cv_class_stratify, fold_id)
  cv_groups_6 <- create_cv_groups(data, dist_6, params, cv_folds, cv_class_stratify, fold_id)
  cv_groups_7 <- create_cv_groups(data, dist_7, params, cv_folds, cv_class_stratify, fold_id)
  cv_groups_8 <- create_cv_groups(data, dist_8, params, cv_folds, cv_class_stratify, fold_id)
  cv_groups_9 <- create_cv_groups(data, dist_9, params, cv_folds, cv_class_stratify, fold_id)
  cv_groups_10 <- create_cv_groups(data, dist_10, params, cv_folds, cv_class_stratify, fold_id)

    
  # Then cv_groups == fold_id
  expect_equal(cv_groups_1, fold_id)
  expect_equal(cv_groups_2, fold_id)
  expect_equal(cv_groups_3, fold_id)
  expect_equal(cv_groups_4, fold_id)
  expect_equal(cv_groups_5, fold_id)
  expect_equal(cv_groups_6, fold_id)
  expect_equal(cv_groups_7, fold_id)
  expect_equal(cv_groups_8, fold_id)
  expect_equal(cv_groups_9, fold_id)
  expect_equal(cv_groups_10, fold_id)
})
test_that("cv groups is equal to fold_id if not NULL - Bernoulli with cv_class_stratify=FALSE", {
  # Given inputs - default distributions for Bernoulli
  # Dist_Objs
  dist <- gbm_dist("Bernoulli")
  
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
  
  # Params Obj
  params <- training_params(num_train = N)
  
  # Other CV Parameters
  cv_folds <- 2
  fold_id <- seq_len(params$num_train)
  
  # When creating cv_groups with cv_class_stratify=FALSE and fold_id != NULL
  cv_class_stratify <- FALSE
  cv_groups <- create_cv_groups(data, dist, params, cv_folds, cv_class_stratify, fold_id)
  
  # Then cv_groups == fold_id
  expect_equal(cv_groups, fold_id)
})
test_that("cv_class_stratify does nothing for distributions other than Bernoulli - fold_id != NULL", {
  # Given inputs - default distributions other than Bernoulli
  # Dist_Objs
  dist_1 <- gbm_dist("AdaBoost")
  dist_2 <- gbm_dist("CoxPH")
  dist_3 <- gbm_dist("Gamma")
  dist_4 <- gbm_dist("Gaussian")
  dist_5 <- gbm_dist("Huberized")
  dist_6 <- gbm_dist("Laplace")
  dist_7 <- gbm_dist("Pairwise")
  dist_8 <- gbm_dist("Poisson")
  dist_9 <- gbm_dist("Quantile")
  dist_10 <- gbm_dist("TDist")
  dist_11 <- gbm_dist("Tweedie")
  
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
  
  # Params Obj
  params <- training_params(num_train = N)
  
  # Other CV Parameters - NB: fold_id not NULL
  cv_folds <- 2
  fold_id <- seq_len(params$num_train)
  
  # When creating cv_groups for both cv_class_stratify==TRUE and FALSE for distributions other than Bernoulli
  cv_groups_1_T <- create_cv_groups(data, dist_1, params, cv_folds, cv_class_stratify=TRUE, fold_id)
  cv_groups_2_T <- create_cv_groups(data, dist_2, params, cv_folds, cv_class_stratify=TRUE, fold_id)
  cv_groups_3_T <- create_cv_groups(data, dist_3, params, cv_folds, cv_class_stratify=TRUE, fold_id)
  cv_groups_4_T <- create_cv_groups(data, dist_4, params, cv_folds, cv_class_stratify=TRUE, fold_id)
  cv_groups_5_T <- create_cv_groups(data, dist_5, params, cv_folds, cv_class_stratify=TRUE, fold_id)
  cv_groups_6_T <- create_cv_groups(data, dist_6, params, cv_folds, cv_class_stratify=TRUE, fold_id)
  cv_groups_7_T <- create_cv_groups(data, dist_7, params, cv_folds, cv_class_stratify=TRUE, fold_id)
  cv_groups_8_T <- create_cv_groups(data, dist_8, params, cv_folds, cv_class_stratify=TRUE, fold_id)
  cv_groups_9_T <- create_cv_groups(data, dist_9, params, cv_folds, cv_class_stratify=TRUE, fold_id)
  cv_groups_10_T <- create_cv_groups(data, dist_10, params, cv_folds, cv_class_stratify=TRUE, fold_id)
  cv_groups_11_T <- create_cv_groups(data, dist_11, params, cv_folds, cv_class_stratify=TRUE, fold_id)
  
  cv_groups_1_F <- create_cv_groups(data, dist_1, params, cv_folds, cv_class_stratify=FALSE, fold_id)
  cv_groups_2_F <- create_cv_groups(data, dist_2, params, cv_folds, cv_class_stratify=FALSE, fold_id)
  cv_groups_3_F <- create_cv_groups(data, dist_3, params, cv_folds, cv_class_stratify=FALSE, fold_id)
  cv_groups_4_F <- create_cv_groups(data, dist_4, params, cv_folds, cv_class_stratify=FALSE, fold_id)
  cv_groups_5_F <- create_cv_groups(data, dist_5, params, cv_folds, cv_class_stratify=FALSE, fold_id)
  cv_groups_6_F <- create_cv_groups(data, dist_6, params, cv_folds, cv_class_stratify=FALSE, fold_id)
  cv_groups_7_F <- create_cv_groups(data, dist_7, params, cv_folds, cv_class_stratify=FALSE, fold_id)
  cv_groups_8_F <- create_cv_groups(data, dist_8, params, cv_folds, cv_class_stratify=FALSE, fold_id)
  cv_groups_9_F <- create_cv_groups(data, dist_9, params, cv_folds, cv_class_stratify=FALSE, fold_id)
  cv_groups_10_F <- create_cv_groups(data, dist_10, params, cv_folds, cv_class_stratify=FALSE, fold_id)
  cv_groups_11_F <- create_cv_groups(data, dist_11, params, cv_folds, cv_class_stratify=FALSE, fold_id)
  
  # These groups don't depend on cv_class_stratify
  expect_equal(cv_groups_1_T, cv_groups_1_F)
  expect_equal(cv_groups_2_T, cv_groups_2_F)
  expect_equal(cv_groups_3_T, cv_groups_3_F)
  expect_equal(cv_groups_4_T, cv_groups_4_F)
  expect_equal(cv_groups_5_T, cv_groups_5_F)
  expect_equal(cv_groups_6_T, cv_groups_6_F)
  expect_equal(cv_groups_7_T, cv_groups_7_F)
  expect_equal(cv_groups_8_T, cv_groups_8_F)
  expect_equal(cv_groups_9_T, cv_groups_9_F)
  expect_equal(cv_groups_10_T, cv_groups_10_F)
  expect_equal(cv_groups_11_T, cv_groups_11_F)
})
test_that("cv_class_stratify does nothing for distributions other than Bernoulli - fold_id == NULL", {
  # Given inputs - default distributions other than Bernoulli
  # Dist_Objs
  dist_1 <- gbm_dist("AdaBoost")
  dist_2 <- gbm_dist("CoxPH")
  dist_3 <- gbm_dist("Gamma")
  dist_4 <- gbm_dist("Gaussian")
  dist_5 <- gbm_dist("Huberized")
  dist_6 <- gbm_dist("Laplace")
  dist_7 <- gbm_dist("Pairwise")
  dist_8 <- gbm_dist("Poisson")
  dist_9 <- gbm_dist("Quantile")
  dist_10 <- gbm_dist("TDist")
  dist_11 <- gbm_dist("Tweedie")
  
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
  
  # Params Obj
  params <- training_params(num_train = N)
  
  # Other CV Parameters - NB: fold_id is NULL
  cv_folds <- 2
  fold_id <- NULL
  
  # When creating cv_groups for both cv_class_stratify==TRUE and FALSE for distributions other than Bernoulli
  set.seed(1)
  cv_groups_1_T <- create_cv_groups(data, dist_1, params, cv_folds, cv_class_stratify=TRUE, fold_id)
  cv_groups_2_T <- create_cv_groups(data, dist_2, params, cv_folds, cv_class_stratify=TRUE, fold_id)
  cv_groups_3_T <- create_cv_groups(data, dist_3, params, cv_folds, cv_class_stratify=TRUE, fold_id)
  cv_groups_4_T <- create_cv_groups(data, dist_4, params, cv_folds, cv_class_stratify=TRUE, fold_id)
  cv_groups_5_T <- create_cv_groups(data, dist_5, params, cv_folds, cv_class_stratify=TRUE, fold_id)
  cv_groups_6_T <- create_cv_groups(data, dist_6, params, cv_folds, cv_class_stratify=TRUE, fold_id)
  cv_groups_7_T <- create_cv_groups(data, dist_7, params, cv_folds, cv_class_stratify=TRUE, fold_id)
  cv_groups_8_T <- create_cv_groups(data, dist_8, params, cv_folds, cv_class_stratify=TRUE, fold_id)
  cv_groups_9_T <- create_cv_groups(data, dist_9, params, cv_folds, cv_class_stratify=TRUE, fold_id)
  cv_groups_10_T <- create_cv_groups(data, dist_10, params, cv_folds, cv_class_stratify=TRUE, fold_id)
  cv_groups_11_T <- create_cv_groups(data, dist_11, params, cv_folds, cv_class_stratify=TRUE, fold_id)
  
  set.seed(1)
  cv_groups_1_F <- create_cv_groups(data, dist_1, params, cv_folds, cv_class_stratify=FALSE, fold_id)
  cv_groups_2_F <- create_cv_groups(data, dist_2, params, cv_folds, cv_class_stratify=FALSE, fold_id)
  cv_groups_3_F <- create_cv_groups(data, dist_3, params, cv_folds, cv_class_stratify=FALSE, fold_id)
  cv_groups_4_F <- create_cv_groups(data, dist_4, params, cv_folds, cv_class_stratify=FALSE, fold_id)
  cv_groups_5_F <- create_cv_groups(data, dist_5, params, cv_folds, cv_class_stratify=FALSE, fold_id)
  cv_groups_6_F <- create_cv_groups(data, dist_6, params, cv_folds, cv_class_stratify=FALSE, fold_id)
  cv_groups_7_F <- create_cv_groups(data, dist_7, params, cv_folds, cv_class_stratify=FALSE, fold_id)
  cv_groups_8_F <- create_cv_groups(data, dist_8, params, cv_folds, cv_class_stratify=FALSE, fold_id)
  cv_groups_9_F <- create_cv_groups(data, dist_9, params, cv_folds, cv_class_stratify=FALSE, fold_id)
  cv_groups_10_F <- create_cv_groups(data, dist_10, params, cv_folds, cv_class_stratify=FALSE, fold_id)
  cv_groups_11_F <- create_cv_groups(data, dist_11, params, cv_folds, cv_class_stratify=FALSE, fold_id)
  
  
  # These groups don't depend on cv_class_stratify
  expect_equal(cv_groups_1_T, cv_groups_1_F)
  expect_equal(cv_groups_2_T, cv_groups_2_F)
  expect_equal(cv_groups_3_T, cv_groups_3_F)
  expect_equal(cv_groups_4_T, cv_groups_4_F)
  expect_equal(cv_groups_5_T, cv_groups_5_F)
  expect_equal(cv_groups_6_T, cv_groups_6_F)
  expect_equal(cv_groups_7_T, cv_groups_7_F)
  expect_equal(cv_groups_8_T, cv_groups_8_F)
  expect_equal(cv_groups_9_T, cv_groups_9_F)
  expect_equal(cv_groups_10_T, cv_groups_10_F)
  expect_equal(cv_groups_11_T, cv_groups_11_F)
  
})
test_that("cv_class_stratify changes grouping for Bernoulli distribution", {
  # Given inputs - default distributions for Bernoulli, data
  # now representative of distribution
  # Dist_Objs
  dist <- gbm_dist("Bernoulli")
  
  # create some data - DOES NEED TO REFLECT DIST
  set.seed(1)
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
  
  offset <- rep(0, N)
  
  X <- data.frame(X1=X1,X2=X2,X3=X3)
  
  data <- gbm_data(X, Y, w, offset)
  
  # Params Obj
  params <- training_params(num_train = N)
  
  # Other CV Parameters
  cv_folds <- 2
  fold_id <- seq_len(params$num_train)
  
  # When creating cv_groups with cv_class_stratify=TRUE and FALSE
  cv_class_stratify_T <- TRUE
  cv_class_stratify_F <- FALSE
  
  set.seed(1)
  cv_groups_T <- create_cv_groups(data, dist, params, cv_folds, cv_class_stratify_T, fold_id)
  
  set.seed(1)
  cv_groups_F <- create_cv_groups(data, dist, params, cv_folds, cv_class_stratify_F, fold_id)
  
  # Then cv_groups are not identical
  expect_true(any(cv_groups_T != cv_groups_F))
})
test_that("Error thrown when number in a certain class is < cv_folds on stratification - Bernoulli", {
  # Given inputs - default distributions for Bernoulli, data
  # now representative of distribution
  # Dist_Objs
  dist <- gbm_dist("Bernoulli")
  
  # create some data - DOES NEED TO REFLECT DIST
  set.seed(1)
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
  
  offset <- rep(0, N)
  
  X <- data.frame(X1=X1,X2=X2,X3=X3)
  
  data <- gbm_data(X, Y, w, offset)
  
  # Params Obj
  params <- training_params(num_train = N)
  
  # Other CV Parameters
  fold_id <- seq_len(params$num_train)
  
  # When creating cv_groups with cv_class_stratify=TRUE and cv_folds > min number of a certain class
  cv_class_stratify <- TRUE
  Ones <- tabulate(data$y[seq_len(params$num_train)])
  Zeros <- length(data$y[seq_len(params$num_train)])-Ones
  cv_folds <- min(c(Ones, Zeros)) + 1
  
  # Then error will be thrown
  expect_error(create_cv_groups(data, dist, params, cv_folds, cv_class_stratify, fold_id))
})
test_that("cv groups are such that folds are split at group boundaries - Pairwise", {
  # Given a Pairwise distribution object and appropriate
  # groupings etc.
  dist <- gbm_dist("Pairwise")

  
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
  
  # Params Obj
  params <- training_params(num_train = N)
  
  # Other CV Parameters - NB: fold_id not NULL
  cv_folds <- 2
  fold_id <- seq_len(params$num_train)
  
  # DEFINE GROUPINGS
  dist$group <- as.factor(sample(seq_len(5), N, replace=TRUE))
  
  # When cv_groups is called and compared to 
  set.seed(1)
  cv_groups <- create_cv_groups(data, dist, params, cv_folds, FALSE, fold_id)
  
  # Then the splits are at group boundaries
  set.seed(1)
  samp <- sample(rep(seq_len(cv_folds), length=nlevels(dist$group)))
  expect_equal(cv_groups, samp[as.integer(dist$group[seq_len(params$num_train)])])
})