####################
# Author: James Hickey
#
# Series of tests to check if plotting functions run correctly
#
####################

context("Testing performance plotting")
test_that("perf_plot runs with all perf methods", {
  # Given a fit object and perf evaluated with all methods
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
  
  fit <- gbm2(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE)
  
  best_iter_t <- best_iter_test(fit)
  best_iter_c <- best_iter_cv(fit)
  best_iter_oo <- best_iter_out_of_bag(fit)
  
  # When calling perf-plot
  # Using all 3 methods
  # Then no error is thrown
  expect_error(perf_plot(fit, best_iter_t, FALSE, FALSE, "test", "main"), NA)
  expect_error(perf_plot(fit, best_iter_t, FALSE, TRUE, "test", "main"), NA)
  expect_error(perf_plot(fit, best_iter_t, TRUE, FALSE, "test", "main"), NA)
  expect_error(perf_plot(fit, best_iter_t, TRUE, TRUE, "test", "main"), NA)
  
  expect_error(perf_plot(fit, best_iter_c, FALSE, FALSE, "cv", "main"), NA)
  expect_error(perf_plot(fit, best_iter_c, FALSE, TRUE, "cv", "main"), NA)
  expect_error(perf_plot(fit, best_iter_c, TRUE, FALSE, "cv", "main"), NA)
  expect_error(perf_plot(fit, best_iter_c, TRUE, TRUE, "cv", "main"), NA)
    
  expect_error(perf_plot(fit, best_iter_oo, FALSE, FALSE, "OOB", "main"), NA)
  expect_error(perf_plot(fit, best_iter_oo, FALSE, TRUE, "OOB", "main"), NA)
  expect_error(perf_plot(fit, best_iter_oo, TRUE, FALSE, "OOB", "main"), NA)
  expect_error(perf_plot(fit, best_iter_oo, TRUE, TRUE, "OOB", "main"), NA)
})

test_that("perf_plot throws error if gbm_fit_obj is not of class GBMFit", {
  # Given a fit object and perf evaluated with all methods
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
  
  fit <- gbm2(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE)
  
  best_iter_t <- best_iter_test(fit)

  # When calling perf_plot with fit stripped of class
  class(fit) <- "Wrong"
  
  # Then no error is thrown
  expect_error(perf_plot(fit, best_iter_t, FALSE, FALSE, "test", "main"))
})

test_that("perf_plot throws error if out_of_bag_curve is not logical or is na", {
  # Given a fit object and perf evaluated with cv method
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
  
  fit <- gbm2(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE)
  
  best_iter_c <- best_iter_cv(fit)
  
  # When calling perf-plot with out_of_bag_curve NA or not logical
  # Then an error is thrown
  expect_error(perf_plot(fit, best_iter_c, out_of_bag_curve=NA, FALSE, "test", "main"))
  expect_error(perf_plot(fit, best_iter_c, out_of_bag_curve=c(TRUE, FALSE), FALSE, "test", "main"))
  expect_error(perf_plot(fit, best_iter_c, out_of_bag_curve=1.5, FALSE, "test", "main"))
  expect_error(perf_plot(fit, best_iter_c, out_of_bag_curve="", FALSE, "test", "main"))
  expect_error(perf_plot(fit, best_iter_c, out_of_bag_curve=NaN, FALSE, "test", "main"))
  
})

test_that("perf_plot throws error if overlay is not logical or is na", {
  # Given a fit object and perf evaluated with cv method
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
  
  fit <- gbm2(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE)
  
  best_iter_c <- best_iter_cv(fit)
  
  # When calling perf-plot with out_of_bag_curve NA or not logical
  # Then an error is thrown
  expect_error(perf_plot(fit, best_iter_c, FALSE, overlay=NA, "test", "main"))
  expect_error(perf_plot(fit, best_iter_c, FALSE, overlay=c(TRUE, FALSE), "test", "main"))
  expect_error(perf_plot(fit, best_iter_c, FALSE, overlay=1.5, "test", "main"))
  expect_error(perf_plot(fit, best_iter_c, FALSE, overlay="", "test", "main"))
  expect_error(perf_plot(fit, best_iter_c, FALSE, overlay=NaN, "test", "main"))
})

test_that("get_ylabel returns correct string for distribution objects", {
  # Given distribution objects
  dist_1 <- gbm_dist("AdaBoost")
  dist_2 <- gbm_dist("Bernoulli")
  dist_3 <- gbm_dist("CoxPH")
  dist_4 <- gbm_dist("Gamma")
  dist_5 <- gbm_dist("Gaussian")
  dist_6 <- gbm_dist("Huberized")
  dist_7 <- gbm_dist("Laplace")
  dist_8_a <- gbm_dist("Pairwise", metric="conc")
  dist_8_b <- gbm_dist("Pairwise", metric="map")
  dist_8_c <- gbm_dist("Pairwise", metric="mrr")
  dist_8_d <- gbm_dist("Pairwise", metric="ndcg")
  dist_9 <- gbm_dist("Poisson")
  dist_10 <- gbm_dist("Quantile")
  dist_11 <- gbm_dist("TDist")
  dist_12 <- gbm_dist("Tweedie")
  
  # When we get the appropriate y_labels for perf_plot
  label_1 <- get_ylabel(dist_1)
  label_2 <- get_ylabel(dist_2)
  label_3 <- get_ylabel(dist_3)
  label_4 <- get_ylabel(dist_4)
  label_5 <- get_ylabel(dist_5)
  label_6 <- get_ylabel(dist_6)
  label_7 <- get_ylabel(dist_7)
  label_8_a <- get_ylabel(dist_8_a)
  label_8_b <- get_ylabel(dist_8_b)
  label_8_c <- get_ylabel(dist_8_c)
  label_8_d <- get_ylabel(dist_8_d)
  label_9 <- get_ylabel(dist_9)
  label_10 <- get_ylabel(dist_10)
  label_11 <- get_ylabel(dist_11)
  label_12 <- get_ylabel(dist_12)
  
  # Then get correct labels
  expect_equal(label_1, "AdaBoost exponential bound")
  expect_equal(label_2, "Bernoulli deviance")
  expect_equal(label_3,"Cox partial deviance")
  expect_equal(label_4, "Gamma deviance")
  expect_equal(label_5, "Squared error loss")
  expect_equal(label_6, "Hinged loss")
  expect_equal(label_7, "Absolute loss")
  expect_equal(label_8_a, "Fraction of concordant pairs")
  expect_equal(label_8_b, "Mean average precision")
  expect_equal(label_8_c, "Mean reciprocal rank")
  expect_equal(label_8_d, "Normalized discounted cumulative gain")
  expect_equal(label_9, "Poisson deviance")
  expect_equal(label_10, "Quantile loss")
  expect_equal(label_11, "t-distribution deviance")
  expect_equal(label_12, "Tweedie deviance")
})

test_that("get_ylabel throws error if passed unrecognised GBMDist class", {
  # Given a GBM dist
  dist <- gbm_dist()
  
  # When GBMDist not recognised
  class(dist) <- c("What is this?", "GBMDist")
  
  # Then error thrown with get_ylabel
  expect_error(get_ylabel(dist))
})

context("Testing calibration plot")
test_that("Error thrown if neither the knots or df specified", {
  
})

test_that("Error thrown if df is not a positive integer", {
  
})

test_that("Error thrown if y and p not same length", {
  
})

test_that("Can run with shade.col not NA", {
  
})

test_that("Can run with defaults", {
  
})

test_that("Can run with replace = TRUE", {
  
})

context("Testinng plot methods for GBMFit")

