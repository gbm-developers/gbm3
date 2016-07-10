####################
# Author: James Hickey
#
# Series of test to check that the conversion of factors 
# and that data is data is appropriate for distribution
#
####################

#### Factor Conversion #### 
context("Testing factor conversion")
test_that("Test factor conversion requires GBMData object", {
  # Given correct data
  N <- 1000
  x <- runif(N)
  p <- 0.5
  y <- rbinom(N,1,p)
  w <- rexp(N) 
  offset <- rexp(N)
  
  # When data object has class removed
  data <- gbm_data(matrix(x), y, w, offset)
  attr(data, "class") <- "FAKE"
  
  # Then error thrown when 
  expect_error(convert_factors(data))
})

#### Validate data given distribution ####
context("Testing data validation and conversion")
test_that("Validation fails if not passed gbm_data object", {
  # Given data (removed class) and a default distribution object
  N <- 1000
  x <- runif(N)
  p <- 0.5
  y <- rbinom(N,1,p)
  w <- rexp(N) 
  offset <- rexp(N)
  
  data <- gbm_data(matrix(x), y, w, offset)
  attr(data, "class") <- "FAKE"
  dist <- gbm_dist()
  
  # Then error thrown
  expect_error(validate_gbm_data(data, dist))
})

test_that("Validation fails if not given distribution obj", {
  # Given data but not a distribution (remove class)
  N <- 1000
  x <- runif(N)
  p <- 0.5
  y <- rbinom(N,1,p)
  w <- rexp(N) 
  offset <- rexp(N)
  
  data <- gbm_data(matrix(x), y, w, offset)
  dist <- gbm_dist()
  attr(dist, "class") <- "FAKE"
  
  # Then error thrown
  expect_error(validate_gbm_data(data, dist))
})

test_that("Weights will be normalized to N if not Pairwise distribution", {
  # Given data and a distribution - not pairwise
  N <- 1000
  x <- runif(N)
  p <- 0.5
  y <- rbinom(N,1,p)
  w <- rexp(N) 
  offset <- rexp(N)
  
  data <- gbm_data(matrix(x), y, w, offset)
  dist <- gbm_dist()
  
  # When data is validated
  data <- validate_gbm_data(data, dist)
  
  # Then weights are normalised to N
  expect_equal(w*length(w)/sum(w), data$weights)
})

test_that("Weights will be normalized across GROUP if Pairwise distribution", {
  # create query groups, with an average size of 25 items each
  set.seed(1)
  N <- 1000
  num.queries <- floor(N/25)
  query <- sample(1:num.queries, N, replace=TRUE)
  
  # X1 is a variable determined by query group only
  query.level <- runif(num.queries)
  X1 <- query.level[query]
  
  # X2 varies with each item
  X2 <- runif(N)
  
  # X3 is uncorrelated with target
  X3 <- runif(N)
  
  # The target
  Y <- X1 + X2
  
  # Add some random noise to X2 that is correlated with
  # queries, but uncorrelated with items
  
  X2 <- X2 + scale(runif(num.queries))[query]
  
  # Add some random noise to target
  SNR <- 5 # signal-to-noise ratio
  sigma <- sqrt(var(Y)/SNR)
  Y <- Y + runif(N, 0, sigma)
  Y.norm <- round(Y/(max(Y) + 0.001))
  
  # Make data and dist
  data <- gbm_data(data.frame(X1, X2), )
  dist <- gbm_dist("Pairwise", metric="ndcg", group=query)
  dist_2 <- gbm_dist("Pairwise", metric="conc", group=query)
  dist_3 <- gbm_dist("Pairwise", metric="mrr", group=query)
  dist_4 <- gbm_dist("Pairwise", metric="map", group=query)
  
  # When validated then weights normalized across GROUP
  
  # Then error thrown on checking responses
  expect_error(check_response_values(dist, Y.norm))
  expect_error(check_response_values(dist_2, Y.norm))
  expect_error(check_response_values(dist_3, Y.norm))
  expect_error(check_response_values(dist_4, Y.norm))
  
})

test_that("Offset vector must contain same number of points as response - if not CoxPH", {
  # Given data and a distribution - not CoxPH
  N <- 1000
  x <- runif(N)
  p <- 0.5
  y <- rbinom(N,1,p)
  w <- rexp(N) 
  offset <- rexp(N)
  
  data <- gbm_data(matrix(x), y, w, offset)
  dist <- gbm_dist()
  
  # When offset does not have the same number of points as the response
  data$offset <- data$offset[1:N-1]
  
  # Then error expected on validation
  expect_error(validate_gbm_data(data, dist))
  
})

test_that("Offset must contain 1/2 number of points as response - if CoxPH", {
  # Given data (not valid but irrelevant here) and a distribution - CoxPH
  # Offset does not have 1/2 number of points as response
  N <- 1000
  x <- runif(N)
  p <- 0.5
  y <- rbinom(N,1,p)
  w <- rexp(N) 
  offset <- rexp(N/2 + 1)
  
  data <- gbm_data(matrix(x), y, w, offset)
  dist <- gbm_dist("CoxPH")
  
  # Then error expected on validation
  expect_error(validate_gbm_data(data, dist))
})

test_that("Responses check requires GBMDist object", {
  # Given data and a distribution - not CoxPH
  N <- 1000
  x <- runif(N)
  p <- 0.5
  y <- rbinom(N,1,p)
  w <- rexp(N) 
  offset <- rexp(N)
  
  data <- gbm_data(matrix(x), y, w, offset)
  dist <- gbm_dist()
  
  # When checking responses without GBMDist object
  attr(dist, "class") <- "FAKE"
  
  # Then error thrown on checking response
  expect_error(check_response_values(dist, data$y))
})

test_that("Responses must be either a data-frame, matrix or vector", {
  # Given data and a distribution - not CoxPH
  N <- 1000
  x <- runif(N)
  p <- 0.5
  y <- rbinom(N,1,p)
  w <- rexp(N) 
  offset <- rexp(N)
  
  data <- gbm_data(matrix(x), y, w, offset)
  dist <- gbm_dist()
  
  # When responses are not a vector/matrix or data-frame
  data$y <- list()
  
  # Then error thrown on checking response
  expect_error(check_response_values(dist, data$y))
})

test_that("AdaBoost responses must be in {0, 1}", {
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  dist <- gbm_dist("AdaBoost")

  # When responses not in {0, 1}
  Y[1] <- -1
  
  # Then error thrown on validation
  expect_error(check_response_values(dist, Y))
})

test_that("Bernoulli responses must be in {0, 1}", {
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  dist <- gbm_dist("Bernoulli")
  
  # When responses not in {0, 1}
  Y[1] <- -1
  
  # Then error thrown on validation
  expect_error(check_response_values(dist, Y))
})

test_that("CoxPH responses must be a survival object", {
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
  dist <- gbm_dist("CoxPH")
  
  # When response is not Surv object
  attr(Resp, "class") <- "NOTSurv"
  
  # Then error thrown on checking response
  expect_error(check_response_values(dist, Resp))
})

test_that("Gamma responses must be positive", {
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  dist <- gbm_dist("Gamma")
  
  # When responses not positive
  Y[1] <- -1
  
  # Then error thrown on validation
  expect_error(check_response_values(dist, Y))
  
})

test_that("Huberized hinge loss requires responses in {0, 1}", {
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  dist <- gbm_dist("Huberized")
  
  # When responses not in {0, 1}
  Y[1] <- -1
  
  # Then error thrown on validation
  expect_error(check_response_values(dist, Y))
  
})

test_that("Pairwise requires non-negative response - all metrics", {
  # create query groups, with an average size of 25 items each
  set.seed(1)
  N <- 1000
  num.queries <- floor(N/25)
  query <- sample(1:num.queries, N, replace=TRUE)
  
  # X1 is a variable determined by query group only
  query.level <- runif(num.queries)
  X1 <- query.level[query]
  
  # X2 varies with each item
  X2 <- runif(N)
  
  # X3 is uncorrelated with target
  X3 <- runif(N)
  
  # The target
  Y <- X1 + X2
  
  # Add some random noise to X2 that is correlated with
  # queries, but uncorrelated with items
  
  X2 <- X2 + scale(runif(num.queries))[query]
  
  # Add some random noise to target
  SNR <- 5 # signal-to-noise ratio
  sigma <- sqrt(var(Y)/SNR)
  Y <- Y + runif(N, 0, sigma)
  Y.norm <- round(Y/(max(Y) + 0.001))
  
  dist <- gbm_dist("Pairwise", metric="ndcg")
  dist_2 <- gbm_dist("Pairwise", metric="conc")
  dist_3 <- gbm_dist("Pairwise", metric="mrr")
  dist_4 <- gbm_dist("Pairwise", metric="map")
  
  # When response is negative
  Y.norm[1] <- -0.01
  
  # Then error thrown on checking responses
  expect_error(check_response_values(dist, Y.norm))
  expect_error(check_response_values(dist_2, Y.norm))
  expect_error(check_response_values(dist_3, Y.norm))
  expect_error(check_response_values(dist_4, Y.norm))
})

test_that("Pairwise map and mrr metrics require response in {0, 1}", {
  # create query groups, with an average size of 25 items each
  set.seed(1)
  N <- 1000
  num.queries <- floor(N/25)
  query <- sample(1:num.queries, N, replace=TRUE)
  
  # X1 is a variable determined by query group only
  query.level <- runif(num.queries)
  X1 <- query.level[query]
  
  # X2 varies with each item
  X2 <- runif(N)
  
  # X3 is uncorrelated with target
  X3 <- runif(N)
  
  # The target
  Y <- X1 + X2
  
  # Add some random noise to X2 that is correlated with
  # queries, but uncorrelated with items
  
  X2 <- X2 + scale(runif(num.queries))[query]
  
  # Add some random noise to target
  SNR <- 5 # signal-to-noise ratio
  sigma <- sqrt(var(Y)/SNR)
  Y <- Y + runif(N, 0, sigma)
  Y.norm <- round(Y/(max(Y) + 0.001))
  
  dist <- gbm_dist("Pairwise", metric="mrr")
  dist_2 <- gbm_dist("Pairwise", metric="map")
  
  # When response is not in {0, 1}
  Y.norm[1] <- 2
  
  # Then error thrown on checking responses
  expect_error(check_response_values(dist, Y.norm))
  expect_error(check_response_values(dist_2, Y.norm))
})

test_that("Poisson requires positive integer response", {
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rpois(N,p)
  dist <- gbm_dist("Poisson")
  
  # When responses not positive integer
  Y[1] <- 0.2
  
  # Then error thrown on validation
  expect_error(check_response_values(dist, Y))
  
})

test_that("Tweedie requires response to be positive", {
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  dist <- gbm_dist("Tweedie")
  
  # When responses not positive
  Y[1] <- -1
  
  # Then error thrown on validation
  expect_error(check_response_values(dist, Y))
})
