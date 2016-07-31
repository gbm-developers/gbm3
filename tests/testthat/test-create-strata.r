####################
# Author: James Hickey
#
# Series of tests to check the functionality that creates strata
#
####################

context("Testing Strata creation:")
test_that("Strata creation function requires GBMData and GBMDist objects", {
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
  offset <- rep(0, N)
  
  Resp <- Surv(tt, delta)
  data <- gbm_data(data.frame(X1, X2, X3), Resp, w, offset)
  train_p <- training_params(id=c(rep(1, N/4), rep(2, N/4), rep(3, N/4), rep(4, N/4)), num_train = 4, 
                             num_features = 3, bag_fraction = 1, min_num_obs_in_node = 1)
  dist <- gbm_dist("CoxPH")
  
  # When not a GBMData object or GBMDist
  copy_data <- data
  copy_dist <- dist
  attr(data, "class") <- "NOTData"
  attr(dist, "class") <- "NOTDist"
  
  # Then an error is thrown
  expect_error(create_strata(data, train_p, copy_dist))
  expect_error(create_strata(copy_data, train_p, dist))
})
test_that("Strata are NA if distribution is not CoxPH", {
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
  offset <- rep(0, N)
  
  Resp <- Surv(tt, delta)
  data <- gbm_data(data.frame(X1, X2, X3), Resp, w, offset)
  train_p <- training_params(id=c(rep(1, N/4), rep(2, N/4), rep(3, N/4), rep(4, N/4)), num_train = 4, 
                             num_features = 3, bag_fraction = 1, min_num_obs_in_node = 1)
  
  # GIVEN Dist Not COXPH
  dist <- gbm_dist("AdaBoost")
  
  # When strata created
  dist <- create_strata(data, train_p, dist)
  
  # Then dist object is unchanged
  expect_true(is.na(dist$strata))
  expect_true(is.na(dist$sorted))
})
test_that("Creating strata fills strata, time_order and sorted fields - CoxPH", {
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
  offset <- rep(0, N/2)
  
  Resp <- Surv(tt, delta)
  data <- gbm_data(data.frame(X1, X2, X3), Resp, w, offset)
  train_p <- training_params(id=c(rep(1, N/4), rep(2, N/4), rep(3, N/4), rep(4, N/4)), num_train = 4, num_features = 3,
                             bag_fraction = 1, min_num_obs_in_node = 1)
  
  # GIVEN Dist - COXPH
  dist <- gbm_dist("CoxPH")
  
  # When strata created
  dist <- create_strata(data, train_p, dist)
  
  # Then dist object is changed
  expect_equal(length(dist$strata), N)
  expect_equal(length(dist$time_order), N)
  expect_equal(nrow(dist$sorted), N)
})
test_that("If response is a matrix with more than 3 columns strata cannot be created - CoxPH", {
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
  offset <- rep(0, N)
  
  
  # GIVEN Dist - COXPH
  dist <- gbm_dist("CoxPH")
  
  # When response has too many columns - set to 4
  Resp <- data.frame(tt, delta)
  Resp <- cbind(cbind(Resp, rnorm(N)), rnorm(N))
  data <- gbm_data(data.frame(X1, X2, X3), Resp, w, offset)
  train_p <- training_params(id=c(rep(1, N/4), rep(2, N/4), rep(3, N/4), rep(4, N/4)), num_train = 4, num_features = 3,
                             bag_fraction = 1, min_num_obs_in_node = 1)
  
  # Then error thrown when creating strata
  expect_error(create_strata(data, train_p, dist))
})
test_that("If strata field in distribution object is NULL, all data are put in same strata - CoxPH", {
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
  offset <- rep(0, N)
  
  Resp <- Surv(tt, delta)
  data <- gbm_data(data.frame(X1, X2, X3), Resp, w, offset)
  train_p <- training_params(id=c(rep(1, N/4), rep(2, N/4), rep(3, N/4), rep(4, N/4)), num_train = 4, num_features = 3, 
                             bag_fraction = 1, min_num_obs_in_node = 1)
  
  # GIVEN Dist - COXPH
  dist <- gbm_dist("CoxPH")
  
  # When creating strata with strata initialized to NULL
  dist$strata <- NULL
  dist <- create_strata(data, train_p, dist)
  
  # Then all examples put in same strata
  expect_equal(dist$strata[1], N)
})
test_that("The training responses are sorted according to strata and this order is stored in time_order", {
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
  offset <- rep(0, N)
  
  Resp <- Surv(tt, delta)
  data <- gbm_data(data.frame(X1, X2, X3), Resp, w, offset)
  train_p <- training_params(id=c(rep(1, N/5), rep(2, N/5), rep(3, N/5), rep(4, N/5)), num_train = 4, num_features = 3, 
                             bag_fraction = 1, min_num_obs_in_node = 1)
  
  # GIVEN Dist - COXPH
  dist <- gbm_dist("CoxPH")
  
  # When strata are created 
  dist <- create_strata(data, train_p, dist)
    
  # Then the training responses are sorted according to strata
  expect_equal(order(-data$y[seq_len(4*N/5), 1]), dist$time_order[seq_len(4*N/5)])
  expect_equal(order(-data$y[(4*N/5 + 1):N, 1])+(4*N/5), dist$time_order[(4*N/5 + 1):N])
})
test_that("Strata not NULL then observations put in different strata - CoxPH", {
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
  offset <- rep(0, N/2)
  
  Resp <- Surv(tt, delta)
  data <- gbm_data(data.frame(X1, X2, X3), Resp, w, offset)
  train_p <- training_params(id=c(rep(1, N/4), rep(2, N/4), rep(3, N/4), rep(4, N/4)), num_train = 4, num_features = 3,
                             bag_fraction = 1, min_num_obs_in_node = 1)
  
  # GIVEN Dist - COXPH
  dist <- gbm_dist("CoxPH")
  Num_Strata <- 5
  strata <- sample(seq_len(Num_Strata), N, replace=TRUE)
  dist$original_strata_id <- strata
  
  # When strata are created
  dist <- create_strata(data, train_p, dist)
  
  # Then ordered by id
  expect_equal(dist$strata[seq_len(Num_Strata)], as.vector(cumsum(table(strata))))
})

