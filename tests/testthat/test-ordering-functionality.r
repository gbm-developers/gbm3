####################
# Author: James Hickey
#
# Series of test to check the data ordering functionality
#
####################

context("Test data ordering functionality")
test_that("Reordering method does nothing to fit unless CoxPH or Pairwise selected", {
  # Given mock gbm_fit
  N <- 1500
  gbm_fit <- structure(list("fit"=1:N), class="GBMFit")
  
  # When ordering according to any distribution other than CoxPH or Pairwise
  dist <- gbm_dist("AdaBoost")
  dist_2 <- gbm_dist("Bernoulli")
  dist_3 <- gbm_dist("Gamma")
  dist_4 <- gbm_dist("Gaussian")
  dist_5 <- gbm_dist("Huberized")
  dist_6 <- gbm_dist("Laplace")
  dist_7 <- gbm_dist("Poisson")
  dist_8 <- gbm_dist("Quantile")
  dist_9 <- gbm_dist("TDist")
  dist_10 <- gbm_dist("Tweedie")
  
  # Then fit order does not change
  expect_equal(reorder_fit(gbm_fit, dist), gbm_fit)
  expect_equal(reorder_fit(gbm_fit, dist_2), gbm_fit)
  expect_equal(reorder_fit(gbm_fit, dist_3), gbm_fit)
  expect_equal(reorder_fit(gbm_fit, dist_4), gbm_fit)
  expect_equal(reorder_fit(gbm_fit, dist_5), gbm_fit)
  expect_equal(reorder_fit(gbm_fit, dist_6), gbm_fit)
  expect_equal(reorder_fit(gbm_fit, dist_7), gbm_fit)
  expect_equal(reorder_fit(gbm_fit, dist_8), gbm_fit)
  expect_equal(reorder_fit(gbm_fit, dist_9), gbm_fit)
  expect_equal(reorder_fit(gbm_fit, dist_10), gbm_fit)
})
test_that("CoxPH- reorders fit according to time and strata", {
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
  Resp <- Surv(tt, delta)
  data <- gbm_data(data.frame(X1, X2, X3), Resp, w, offset=rep(0, N))
  train_p <- training_params(id=c(rep(1, N/5), rep(2, N/5), rep(3, N/5), rep(4, N/5), rep(5, N/5)), num_train = 5, num_features = 3,
                             bag_fraction=1, min_num_obs_in_node=1)
  dist <- gbm_dist("CoxPH")
  
  # Reorder Fit
  dist <- create_strata(data, train_p, dist)
  gbm_fit <- structure(list("fit"=1:N), class="GBMFit")
  gbm_fit_ordered <- reorder_fit(gbm_fit, dist)
  
  # Expect order according to time_order
  gbm_fit$fit[dist$time_order] <- gbm_fit$fit
  expect_equal(gbm_fit_ordered$fit, gbm_fit$fit)
})
test_that("Pairwise - reorders fit according to group order", {
  # Given Pairwise data and fit
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
  dist$group_order <-order(dist$group, -Y)
  fit <- structure(list("fit"=1:N), class="GBMFit")
  
  # When fit is reordered
  reordered_fit <- reorder_fit(fit, dist)
  
  # Then correctly reorders fit
  expect_equal(reordered_fit$fit, fit$fit[order(dist$group_order)])
})
test_that("Order data requires a gbm_data object and distribution object", {
  # Given gbm_data and dist
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
  train_p <- training_params(id=rep(1, N), num_train = N, num_features = 3)
  dist <- gbm_dist("CoxPH")
  
  # When classes is removed
  data_2 <- data
  dist_2 <- dist
  attr(data, "class") <- "NOTData"
  attr(dist, "class") <- "NOTDist"
  
  # Then error thrown when trying to order data
  expect_error(order_data(data, train_p, dist_2))
  expect_error(order_data(data_2, train_p, dist))
})
test_that("Can order data by id", {
  # Given GBM data and training parameters
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
  train_p <- training_params(id=sample(1:5, N, replace=TRUE), num_train = 5, bag_fraction=1, min_num_obs_in_node = 1, num_features = 3)
  
  # When ordered by id
  data_copy <- data
  data <- order_by_id(data, train_p)
  
  # Then data is ordered by id
  expect_equal(data$x, data_copy$x[train_p$id, , drop=FALSE])
  expect_equal(data$y, data_copy$y[train_p$id])
})
test_that("Can order data by groupings", {
  # Given pairwise data and dist
  # Given Pairwise data and fit
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
  dist$group_order <-order(dist$group, -Y)
  
  
  # When order_by_groupings is called
  data_ordered <- order_by_groupings(data, dist)
  
  # Then correctly orders by groupings
  expect_equal(data_ordered$x, data$x[dist$group_order,,drop=FALSE])
  expect_equal(data_ordered$y, data$y[dist$group_order])
  expect_equal(data_ordered$weights, data$weights[dist$group_order])
  
})
test_that("Ordering by groupings does nothing if not Pairwise", {
  # Given data and not pairwise (using CoxPH data for example)
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
  train_p <- training_params(id=sample(1:5, N, replace=TRUE), num_train = 5, bag_fraction=1, min_num_obs_in_node = 1,
                             num_features = 3)
  
  # When ordered by grouping - NOT Pairwise dist
  # Then nothing changes
  expect_equal(data, order_by_groupings(data, gbm_dist("AdaBoost")))
  expect_equal(data, order_by_groupings(data, gbm_dist("Bernoulli")))
  expect_equal(data, order_by_groupings(data, gbm_dist("Gamma")))
  expect_equal(data, order_by_groupings(data, gbm_dist("Gaussian")))
  expect_equal(data, order_by_groupings(data, gbm_dist("Huberized")))
  expect_equal(data, order_by_groupings(data, gbm_dist("Laplace")))
  expect_equal(data, order_by_groupings(data, gbm_dist("Poisson")))
  expect_equal(data, order_by_groupings(data, gbm_dist("Quantile")))
  expect_equal(data, order_by_groupings(data, gbm_dist("TDist")))
  expect_equal(data, order_by_groupings(data, gbm_dist("Tweedie")))
})
test_that("Can order by predictors - check predictor_order", {
  # Given gbm data and training parameters
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
  train_p <- training_params(id=sample(1:5, N, replace=TRUE), num_train = 5, bag_fraction = 1, min_num_obs_in_node = 1, num_features = 3)
  
  # Then can order by predictors
  data <- predictor_order(data, train_p)
  expect_equal(data$x_order, apply(data$x[seq_len(train_p$num_train_rows),,drop=FALSE], 2, order,na.last=FALSE)-1)
})