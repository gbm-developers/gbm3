####################
# Author: James Hickey
#
# Series of test to validate the GBMVarCont objects
#
####################


##### Definition #####
context("Testing Variable Properties container definition:")
test_that("Constructed object has correct class", {
  # Given correct gbm-data object
  N <- 1000
  x <- runif(N)
  p <- 0.5
  y <- rbinom(N,1,p)
  w <- rexp(N) 
  offset <- rexp(N)
  data <- gbm_data(matrix(x), y, w, offset)
  
  # Constructed container has correct class
  expect_true("GBMVarCont" %in% class(var_container(data)))
})

test_that("Constructed object has correct fields", {
  # Given correct gbm-data object
  N <- 1000
  x <- runif(N)
  p <- 0.5
  y <- rbinom(N,1,p)
  w <- rexp(N) 
  offset <- rexp(N)
  data <- gbm_data(matrix(x), y, w, offset)
  
  # Constructed container has correct fields
  expect_equal(names(var_container(data)), c("var_monotone", "var_names", "var_levels", "var_type"))
})

test_that("If variable names are not provided they and retrieved from data", {
  # Given correct gbm-data object - variables have names
  N <- 1000
  x <- runif(N)
  p <- 0.5
  y <- rbinom(N,1,p)
  w <- rexp(N) 
  offset <- rexp(N)
  
  x <- matrix(x)
  colnames(x) <- "Test"
  data <- gbm_data(x, y, w, offset)
  
  # Then constructed container extracts colnames
  expect_equal(var_container(data)$var_names, "Test")
})

test_that("If monotonicity is not provided it defaults to 0", {
  # Given correct gbm-data object
  N <- 1000
  x <- runif(N)
  p <- 0.5
  y <- rbinom(N,1,p)
  w <- rexp(N) 
  offset <- rexp(N)
  data <- gbm_data(matrix(x), y, w, offset)
  
  # Constructed container will have monotonicity of 0 for
  # each predictor
  expect_equal(var_container(data)$var_monotone, rep(0, ncol(matrix(x))) )
})

##### Errors ####
context("Testing incorrect inputs will throw errors:")
test_that("Names of data must be vector of strings", {
  # Given correct gbm-data object
  N <- 1000
  x <- matrix(c(runif(N), runif(N)), nrow=N, ncol=2)
  p <- 0.5
  y <- rbinom(N,1,p)
  w <- rexp(N) 
  offset <- rexp(N)
  data <- gbm_data(x, y, w, offset)
  
  # When names are not vector of strings
  names <- list("1", "2")
  names_2 <- c(NULL, "2")
  names_3 <- c(NA, "2")
  
  # Then error is thrown on construction
  expect_error(var_container(data, var_names = names))
  expect_error(var_container(data, var_names = names_2))
  expect_error(var_container(data, var_names = names_3))
})

test_that("Monotonicity must be specified for each predictor", {
  # Given correct gbm-data object
  N <- 1000
  x <- matrix(c(runif(N), runif(N)), nrow=N, ncol=2)
  p <- 0.5
  y <- rbinom(N,1,p)
  w <- rexp(N) 
  offset <- rexp(N)
  data <- gbm_data(x, y, w, offset)
  
  # When monotonicity not specified for each predictor
  monotone <- c(-1)
  
  # Then expect error thrown
  expect_error(var_container(data, var_monotone = monotone))
})

test_that("Monotonicity must be -1, 0 or 1", {
  # Given correct gbm-data object
  N <- 1000
  x <- matrix(c(runif(N), runif(N)), nrow=N, ncol=2)
  p <- 0.5
  y <- rbinom(N,1,p)
  w <- rexp(N) 
  offset <- rexp(N)
  data <- gbm_data(x, y, w, offset)
  
  # When monotonicity is not 0, -1, +1
  monotone <- c(-1, NA)
  monotone_2 <- c(0.5, 1)
  monotone_3 <- c(-2, 4)
  monotone_4 <- c(0, Inf)
  
  # Then expect error thrown
  expect_error(var_container(data, var_monotone = monotone))
  expect_error(var_container(data, var_monotone = monotone_2))
  expect_error(var_container(data, var_monotone = monotone_3))
  expect_error(var_container(data, var_monotone = monotone_4))
})

test_that("The number of variable names MUST match the number of predictors", {
  # Given correct gbm-data object
  N <- 1000
  x <- runif(N)
  p <- 0.5
  y <- rbinom(N,1,p)
  w <- rexp(N) 
  offset <- rexp(N)
  data <- gbm_data(matrix(x), y, w, offset)
  
  # When number of var_names does not match the number of predictors
  var_names <- c("First", "OneTooMany")
  
  # Then error on construction
  expect_error(var_container(data, var_names = var_names))
})