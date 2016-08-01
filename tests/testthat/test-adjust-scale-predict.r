####################
# Author: James Hickey
#
# Series of tests to check if prediction scales are adjusted correctly.
#
####################

context("Testing adjust_pred_scale")
test_that("AdaBoost predictions are adjusted correctly", {
  # Given a set of predictions and the appropriate distribution obj
  N <- 100
  preds <- runif(N)
  dist <- gbm_dist("AdaBoost")
  
  # When predictions are adjusted
  adj_preds <- adjust_pred_scale(preds, dist)
  
  # Then they are correctly adjusted
  expect_equal(adj_preds, 1/(1 + exp(-2 * preds)))
})
test_that("Bernoulli predictions are adjusted correctly", {
  # Given a set of predictions and the appropriate distribution obj
  N <- 100
  preds <- runif(N)
  dist <- gbm_dist("Bernoulli")
  
  # When predictions are adjusted
  adj_preds <- adjust_pred_scale(preds, dist)
  
  # Then they are correctly adjusted
  expect_equal(adj_preds, 1/(1 + exp(-preds)))
})
test_that("CoxPH predictions are adjusted correctly", {
  # Given a set of predictions and the appropriate distribution obj
  N <- 100
  preds <- runif(N)
  dist <- gbm_dist("CoxPH")
  
  # When predictions are adjusted
  adj_preds <- adjust_pred_scale(preds, dist)
  
  # Then they are correctly adjusted
  expect_equal(adj_preds, preds)
})
test_that("Gamma predictions are adjusted correctly", {
  # Given a set of predictions and the appropriate distribution obj
  N <- 100
  preds <- runif(N)
  dist <- gbm_dist("Gamma")
  
  # When predictions are adjusted
  adj_preds <- adjust_pred_scale(preds, dist)
  
  # Then they are correctly adjusted
  expect_equal(adj_preds, exp(preds))
})
test_that("Gaussian predictions are adjusted correctly", {
  # Given a set of predictions and the appropriate distribution obj
  N <- 100
  preds <- runif(N)
  dist <- gbm_dist("Gaussian")
  
  # When predictions are adjusted
  adj_preds <- adjust_pred_scale(preds, dist)
  
  # Then they are correctly adjusted
  expect_equal(adj_preds, preds)
})
test_that("Huberized predictions are adjusted correctly", {
  # Given a set of predictions and the appropriate distribution obj
  N <- 100
  preds <- runif(N)
  dist <- gbm_dist("Huberized")
  
  # When predictions are adjusted
  adj_preds <- adjust_pred_scale(preds, dist)
  
  # Then they are correctly adjusted
  expect_equal(adj_preds, preds)
})
test_that("Laplace predictions are adjusted correctly", {
  # Given a set of predictions and the appropriate distribution obj
  N <- 100
  preds <- runif(N)
  dist <- gbm_dist("Laplace")
  
  # When predictions are adjusted
  adj_preds <- adjust_pred_scale(preds, dist)
  
  # Then they are correctly adjusted
  expect_equal(adj_preds, preds)
})
test_that("Pairwise predictions are adjusted correctly", {
  # Given a set of predictions and the appropriate distribution obj
  N <- 100
  preds <- runif(N)
  dist <- gbm_dist("Pairwise")
  
  # When predictions are adjusted
  adj_preds <- adjust_pred_scale(preds, dist)
  
  # Then they are correctly adjusted
  expect_equal(adj_preds, 1/(1 + exp(-preds)))
})
test_that("Poisson predictions are adjusted correctly", {
  # Given a set of predictions and the appropriate distribution obj
  N <- 100
  preds <- runif(N)
  dist <- gbm_dist("Poisson")
  
  # When predictions are adjusted
  adj_preds <- adjust_pred_scale(preds, dist)
  
  # Then they are correctly adjusted
  expect_equal(adj_preds, exp(preds))
})
test_that("Quantile predictions are adjusted correctly", {
  # Given a set of predictions and the appropriate distribution obj
  N <- 100
  preds <- runif(N)
  dist <- gbm_dist("Quantile")
  
  # When predictions are adjusted
  adj_preds <- adjust_pred_scale(preds, dist)
  
  # Then they are correctly adjusted
  expect_equal(adj_preds, preds)
})
test_that("TDist predictions are adjusted correctly", {
  # Given a set of predictions and the appropriate distribution obj
  N <- 100
  preds <- runif(N)
  dist <- gbm_dist("TDist")
  
  # When predictions are adjusted
  adj_preds <- adjust_pred_scale(preds, dist)
  
  # Then they are correctly adjusted
  expect_equal(adj_preds, preds)
})
test_that("Tweedie predictions are adjusted correctly", {
  # Given a set of predictions and the appropriate distribution obj
  N <- 100
  preds <- runif(N)
  dist <- gbm_dist("Tweedie")
  
  # When predictions are adjusted
  adj_preds <- adjust_pred_scale(preds, dist)
  
  # Then they are correctly adjusted
  expect_equal(adj_preds,  exp(preds))
})
