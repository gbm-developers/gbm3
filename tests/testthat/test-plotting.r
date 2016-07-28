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

context("Testing quantile_rug plot")
test_that("quantile_rug runs without error", {
  probs <- rnorm(100)
  probs <- abs(probs)/sum(abs(probs))
  
  expect_error(quantile_rug(rnorm(100)), NA)
  expect_error(quantile_rug(rnorm(100), probs), NA)
})
test_that("quantile_rug throws error if probabilities outsider [0, 1]", {
  probs <- rnorm(100)
  probs[1] <- 5
  expect_error(quantile_rug(rnorm(100), probs))
})
test_that("output from quantile_rug is correct - no jittering", {
  probs <- (0:10)/10
  x <- rnorm(100)
  output <- quantile_rug(x, probs)
  expect_equal(output, rug(quantile(x[!is.na(x)], probs)))
})
test_that("quantile_rug jitters the inputs if quantiles < length(probabilities)", {
  set.seed(1)
  x <- rep(1:5, 2)
  probs <- (0:100)/100
  output <- quantile_rug(x, probs)
  
  set.seed(1)
  true_output <- quantile(x, probs)
  true_output <- jitter(true_output)
  
  expect_equal(output, rug(true_output))
})

context("Testing calibration plot")
test_that("Error thrown if neither the knots or df specified - both NULL", {
  # Given data - from example
  require(rpart)
  data(kyphosis)
  y <- as.numeric(kyphosis$Kyphosis)-1
  x <- kyphosis$Age
  glm1 <- glm(y~poly(x,2),family=binomial)
  p <- predict(glm1, type="response")
  
  # Then error thrown when both knots and df are NULL
  expect_error(calibrate_plot(y, p, df=NULL, knots=NULL, xlim=c(0,0.6), ylim=c(0,0.6)))
})
test_that("Error thrown if df is not a positive integer (if vector first element must be a positive integer)", {
  # Given data - from example
  require(rpart)
  data(kyphosis)
  y <- as.numeric(kyphosis$Kyphosis)-1
  x <- kyphosis$Age
  glm1 <- glm(y~poly(x,2),family=binomial)
  p <- predict(glm1, type="response")
  
  # Then error thrown when df Not a positive integer
  expect_error(calibrate_plot(y, p, df=0, knots=NULL, xlim=c(0,0.6), ylim=c(0,0.6)))
  expect_error(calibrate_plot(y, p, df=1.4, knots=NULL, xlim=c(0,0.6), ylim=c(0,0.6)))
  expect_error(calibrate_plot(y, p, df="Wrong", knots=NULL, xlim=c(0,0.6), ylim=c(0,0.6)))
  expect_error(calibrate_plot(y, p, df=c(1, 2), knots=NULL, xlim=c(0,0.6), ylim=c(0,0.6)), NA)
  expect_error(calibrate_plot(y, p, df=c(1.4, 2), knots=NULL, xlim=c(0,0.6), ylim=c(0,0.6)))
  expect_error(calibrate_plot(y, p, df=NaN, knots=NULL, xlim=c(0,0.6), ylim=c(0,0.6)))
  expect_error(calibrate_plot(y, p, df=NA, knots=NULL, xlim=c(0,0.6), ylim=c(0,0.6)))
  expect_error(calibrate_plot(y, p, df=Inf, knots=NULL, xlim=c(0,0.6), ylim=c(0,0.6)))
})
test_that("Error thrown if y and p not same length", {
  # Given data - from example - but y and p not same length now
  require(rpart)
  data(kyphosis)
  y <- as.numeric(kyphosis$Kyphosis)-1
  x <- kyphosis$Age
  glm1 <- glm(y~poly(x,2),family=binomial)
  p <- predict(glm1, type="response")
  p <- p[seq(length(p)-1)]
  
  # Then error thrown
  expect_error(calibrate_plot(y, p))
})
test_that("Can run with defaults", {
  # Given data - from example - but y and p not same length now
  require(rpart)
  data(kyphosis)
  y <- as.numeric(kyphosis$Kyphosis)-1
  x <- kyphosis$Age
  glm1 <- glm(y~poly(x,2),family=binomial)
  p <- predict(glm1, type="response")
  
  
  # Then no error thrown
  expect_error(calibrate_plot(y, p), NA)
})
test_that("Can run with shade_col not NA", {
  # Given data - from example - but y and p not same length now
  require(rpart)
  data(kyphosis)
  y <- as.numeric(kyphosis$Kyphosis)-1
  x <- kyphosis$Age
  glm1 <- glm(y~poly(x,2),family=binomial)
  p <- predict(glm1, type="response")
  
  
  # Then no error thrown
  expect_error(calibrate_plot(y, p, shade_col=1), NA)
})
test_that("Can run with replace = FALSE", {
  # Given data - from example - but y and p not same length now
  require(rpart)
  data(kyphosis)
  y <- as.numeric(kyphosis$Kyphosis)-1
  x <- kyphosis$Age
  glm1 <- glm(y~poly(x,2),family=binomial)
  p <- predict(glm1, type="response")
  
  
  # Then no error thrown
  expect_error(calibrate_plot(y, p, replace=FALSE), NA)
})
test_that("Can run with shade_density != NULL", {
  # Given data - from example - but y and p not same length now
  require(rpart)
  data(kyphosis)
  y <- as.numeric(kyphosis$Kyphosis)-1
  x <- kyphosis$Age
  glm1 <- glm(y~poly(x,2),family=binomial)
  p <- predict(glm1, type="response")
  
  
  # Then no error thrown
  expect_error(calibrate_plot(y, p, shade_density=2.0), NA)
})
test_that("Can run  all distributions", {
  # Given data - from example - but y and p not same length now
  require(rpart)
  data(kyphosis)
  y <- as.numeric(kyphosis$Kyphosis)-1
  x <- kyphosis$Age
  glm1 <- glm(y~poly(x,2),family=binomial)
  p <- predict(glm1, type="response")
  

  # Then no error thrown
  expect_error(calibrate_plot(y, p, distribution = "AdaBoost"), NA)
  expect_error(calibrate_plot(y, p, distribution = "Bernoulli"), NA)
  expect_error(calibrate_plot(y, p, distribution = "CoxPH"), NA)
  expect_error(calibrate_plot(y, p, distribution = "Gamma"), NA)
  expect_error(calibrate_plot(y, p, distribution = "Gaussian"), NA)
  expect_error(calibrate_plot(y, p, distribution = "Laplace"), NA)
  expect_error(calibrate_plot(y, p, distribution = "Huberized"), NA)
  expect_error(calibrate_plot(y, p, distribution = "Pairwise"), NA)
  expect_error(calibrate_plot(y, p, distribution = "Poisson"), NA)
  expect_error(calibrate_plot(y, p, distribution = "Quantile"), NA)
  expect_error(calibrate_plot(y, p, distribution = "TDist"), NA)
  expect_error(calibrate_plot(y, p, distribution = "Tweedie"), NA)
})

context("Testinng plot methods for GBMFit")
test_that("Error thrown if type of plot not 'link' or 'response' ", {
  
})

test_that("Error thrown if var_index has variable outside range not used in fit", {
  
})

test_that("Response throws warning if not Bernoulli, Poisson, Gamma, Pairwise or Tweedie", {
  # Given distribution and response
  dist_1 <- gbm_dist("AdaBoost")
  dist_2 <- gbm_dist("CoxPH")
  dist_3 <- gbm_dist("Gaussian")
  dist_4 <- gbm_dist("Huberized")
  dist_5 <- gbm_dist("Laplace")
  dist_6 <- gbm_dist("Quantile")
  dist_7 <- gbm_dist("TDist")
  resp <- rep(1, 10)
  
  # When response S3 method called
  # Then warning thrown
  expect_warning(response(resp, dist_1))
  expect_warning(response(resp, dist_2))
  expect_warning(response(resp, dist_3))
  expect_warning(response(resp, dist_4))
  expect_warning(response(resp, dist_5))
  expect_warning(response(resp, dist_6))
  expect_warning(response(resp, dist_7))
  
})
test_that("Response returns NULL if not Bernoulli, Poisson, Gamma, Pairwise or Tweedie", {
  # Given distribution and response
  dist_1 <- gbm_dist("AdaBoost")
  dist_2 <- gbm_dist("CoxPH")
  dist_3 <- gbm_dist("Gaussian")
  dist_4 <- gbm_dist("Huberized")
  dist_5 <- gbm_dist("Laplace")
  dist_6 <- gbm_dist("Quantile")
  dist_7 <- gbm_dist("TDist")
  resp <- rep(1, 10)
  
  # When response S3 method called
  # Then returns NULL
  expect_null(response(resp, dist_1))
  expect_null(response(resp, dist_2))
  expect_null(response(resp, dist_3))
  expect_null(response(resp, dist_4))
  expect_null(response(resp, dist_5))
  expect_null(response(resp, dist_6))
  expect_null(response(resp, dist_7))
})
test_that("Response returns correct answer for Bernoulli, Pairwise, Poisson, Gamma and Tweedie", {
  # Given distributions which don't return NULL and response
  dist_1 <- gbm_dist("Bernoulli")
  dist_2 <- gbm_dist("Pairwise")
  dist_3 <- gbm_dist("Poisson")
  dist_4 <- gbm_dist("Gamma")
  dist_5 <- gbm_dist("Tweedie")
  
  resp <- rep(1, 10)
  
  # When response S3 method called
  # Then returns correct answer
  expect_equal(response(resp, dist_1), 1/(1+exp(-resp)))
  expect_equal(response(resp, dist_2), 1/(1+exp(-resp)))
  expect_equal(response(resp, dist_3), exp(resp))
  expect_equal(response(resp, dist_4), exp(resp))
  expect_equal(response(resp, dist_5), exp(resp))
})

test_that("get_default_grid_levels returns answer of correct type", {
  
})

test_that("generate_grid_levels throws an error if length of grid_levels not same length as var_index", {
  
})

test_that("generate_grid_levels can run without error", {
  
})

test_that("warning thrown if num_trees exceeds those in fit", {
  
})

test_that("warning thrown if number var indices > 3", {
  
})

test_that("return_grid=TRUE returns the grid", {
  
})

test_that("number of var indices >3 sets return_grid to TRUE", {
  
})

test_that("can plot with one variable selected", {
  
})

test_that("can correctly get one variable y-label", {
  
})

test_that("can plot with two variables", {
  
})

test_that("can plot with 3 variables selected", {
  
})