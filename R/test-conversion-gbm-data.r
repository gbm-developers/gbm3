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
  
})


#### Validate data given distribution ####
context("Testing data validation and conversion")
test_that("Validation fails if not passed gbm_data object", {
  # Given data (removed class) and a default distribution object
  
  # Then error thrown
})

test_that("Validation fails if not given distribution obj", {
  # Given data but not a distribution (remove class)
  
  # Then error thrown
})


test_that("Weights must be doubles", {
  
})

test_that("Weights must be positive", {
  
})

test_that("Weights will be normalized if not Pairwise distribution", {
  
})

test_that("Weights will be normalized across GROUP if Pairwise distribution", {
  
})

test_that("Offset vector must contain same number of points as response - if not CoxPH", {
  
})

test_that("Offset must contain 1/2 number of points as response - if CoxPH", {
  
})

test_that("Responses check requires GBMDistObj", {
  
})

test_that("Responses must be either a data-frame, matrix or vector", {
  
})

test_that("AdaBoost responses must be in {0, 1}", {
  
})

test_that("Bernoulli responses must be in {0, 1}", {
  
})

test_that("CoxPH responses must be a survival object", {
  
})

test_that("Gamma responses must be positive", {
  
})

test_that("Huberized hinge loss requires responses in {0, 1}", {
  
})

test_that("Pairwise requires non-negative response - all metrics", {
  
})

test_that("Pairwise map and mrr metrics require response in {0, 1}", {
  
})

test_that("Poisson requires positive integer response", {
  
})

test_that("Tweedie requires response to be positive", {
  
})
