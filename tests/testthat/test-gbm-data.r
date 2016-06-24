####################
# Author: James Hickey
#
# Series of test to validate the GBMData objects
#
####################

##### Definition #####
context("Testing GBMData object definition:")
test_that("Constructed object has correct class", {
  # Given correct data
  N <- 1000
  x <- runif(N)
  p <- 0.5
  y <- rbinom(N,1,p)
  w <- rexp(N) 
  offset <- rexp(N)
  
  # Constructed object
  expect_true("GBMData" %in% class(gbm_data(matrix(x), y, w, offset)))
})

test_that("Constructed object has correct fields", {
  # Given correct data
  N <- 1000
  x <- runif(N)
  p <- 0.5
  y <- rbinom(N,1,p)
  w <- rexp(N) 
  offset <- rexp(N)
  
  # Constructed objetc
  expect_equal(names(gbm_data(matrix(x), y, w, offset)), c("x", "y", "weights", "offset"))
})

#### Error Checking #####
context("Testing variable check")
test_that("Predictors must be in a data-frame or matrix", {
  # Given correct data
  N <- 1000
  x <- runif(N)
  p <- 0.5
  y <- rbinom(N,1,p)
  w <- rexp(N) 
  offset <- rexp(N)
  
  # Error thrown when predictors not put in matrix or data-frame
  expect_error(gbm_data(x, y, w, offset))
  expect_error(gbm_data(matrix(x), y, w, offset), NA)
  expect_error(gbm_data(data.frame(x), y, w, offset), NA)
})

test_that("Responses must be in a matrix or vector", {
  # Given correct data
  N <- 1000
  x <- runif(N)
  p <- 0.5
  y <- rbinom(N,1,p)
  w <- rexp(N) 
  offset <- rexp(N)
  
  # Error thrown when responses not put in matrix or vector
  expect_error(gbm_data(matrix(x), list(y), w, offset))
  expect_error(gbm_data(matrix(x), y, w, offset), NA)
  expect_error(gbm_data(matrix(x), matrix(y), w, offset), NA)
  expect_error(gbm_data(matrix(x), data.frame(y), w, offset))
})

test_that("Each predictor row must have a response", {
  # Given correct data
  N <- 1000
  x <- runif(N-1)
  x2 <- runif(N+1)
  p <- 0.5
  y <- rbinom(N,1,p)
  w <- rexp(N) 
  offset <- rexp(N)
  
  # Error thrown as more responses than predictors and vice-versa
  expect_error(gbm_data(matrix(x), y, w, offset))
  expect_error(gbm_data(matrix(x2), y, w, offset))
})

test_that("Predictor variables must be ordered, numeric or categorical", {
  # Given correct data
  N <- 1000
  x <- runif(N)
  p <- 0.5
  y <- rbinom(N,1,p)
  w <- rexp(N) 
  offset <- rexp(N)
  
  # Error thrown when x is not categorical, ordered or numeric
  expect_error(gbm_data(matrix(rep("A", N)), y, w, offset))
  expect_error(gbm_data(data.frame(ordered(x)), y, w, offset), NA)
  expect_error(gbm_data(data.frame(factor(x)), y, w, offset), NA)
  expect_error(gbm_data(matrix(as.numeric(x)), y, w, offset), NA)
})

test_that("Missing values not allowed in response", {
  # Given correct data
  N <- 1000
  x <- runif(N)
  p <- 0.5
  y <- rbinom(N,1,p)
  w <- rexp(N) 
  offset <- rexp(N)
  
  # Set subset of responses to be Missing
  y[sample(1:N, 100)] <- NA
  
  # Expect error
  expect_error(gbm_data(matrix(x), y, w, offset))
})

test_that("Can't have predictor with all missing values", {
  # Given correct data
  N <- 1000
  x <- runif(N)
  p <- 0.5
  y <- rbinom(N,1,p)
  w <- rexp(N) 
  offset <- rexp(N)
  
  # Add in variable which is all NAs
  x2 <- rep(NA, N)
  
  # Expect error to be thrown
  expect_error(gbm_data(data.frame(x, x2), y, w, offset))
})

test_that("Weights must be a vector of positive doubles", {
  # Given correct data
  N <- 1000
  x <- runif(N)
  p <- 0.5
  y <- rbinom(N,1,p)
  offset <- rexp(N)
  
  expect_error(gbm_data(matrix(x), y, w=-rexp(N), offset))
})

test_that("Offset must be a vector of doubles", {
  # Given correct data
  N <- 1000
  x <- runif(N)
  p <- 0.5
  y <- rbinom(N, 1, p)
  w <- rexp(N)
  
  expect_error(gbm_data(matrix(x), y, w, c(1, "A")))
  expect_error(gbm_data(matrix(x), y, w, 1))
  expect_error(gbm_data(matrix(x), y, w, Inf))
  expect_error(gbm_data(matrix(x), y, w, NULL))
})