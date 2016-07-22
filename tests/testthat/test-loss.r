####################
# Author: James Hickey
#
# Series of tests for gbm loss functionality
#
####################

context("Error checking for GBM loss methods")
test_that("Error thrown responses not vector of doubles", {
  # Given offset, weights, predictions and default dist
  N <- 100
  dist <- gbm_dist()
  offset <- runif(N)
  weights <- runif(N)
  preds <- runif(N)
  
  # When responses are not a vector of doubles
  resps <- rep("1", N)
  resps_2 <- NA
  resps_3 <- Inf
  resps_4 <- rep(NaN, N)
  
  # Then error thrown when trying to calculate the loss
  expect_error(loss(resps, preds, weights, offset, dist))
  expect_error(loss(resps_2, preds, weights, offset, dist))
  expect_error(loss(resps_3, preds, weights, offset, dist))
  expect_error(loss(resps_4, preds, weights, offset, dist))
})
test_that("Error thrown when predictions not vector of doubles", {
  # Given offset, weights, responses and default dist
  N <- 100
  dist <- gbm_dist()
  offset <- runif(N)
  weights <- runif(N)
  resps <- runif(N)
  
  # When predictions are not a vector of doubles
  preds <- rep("1", N)
  preds_2 <- NA
  preds_3 <- Inf
  preds_4 <- rep(NaN, N)
  
  # Then error thrown when trying to calculate the loss
  expect_error(loss(resps, preds, weights, offset, dist))
  expect_error(loss(resps, preds_2, weights, offset, dist))
  expect_error(loss(resps, preds_3, weights, offset, dist))
  expect_error(loss(resps, preds_4, weights, offset, dist))
})
test_that("Error thrown when weights not vector of doubles", {
  # Given offset, responses, predictions and default dist
  N <- 100
  dist <- gbm_dist()
  offset <- runif(N)
  resps <- runif(N)
  preds <- runif(N)
  
  # When weights are not a vector of doubles
  weights <- rep("1", N)
  weights_2 <- NA
  weights_3 <- Inf
  weights_4 <- rep(NaN, N)
  
  # Then error thrown when trying to calculate the loss
  expect_error(loss(resps, preds, weights, offset, dist))
  expect_error(loss(resps, preds, weights_2, offset, dist))
  expect_error(loss(resps, preds, weights_3, offset, dist))
  expect_error(loss(resps, preds, weights_4, offset, dist))
})
test_that("Error thrown when offset not vector of doubles", {
  # Given responses, weights, predictions and default dist
  N <- 100
  dist <- gbm_dist()
  resps <- runif(N)
  weights <- runif(N)
  preds <- runif(N)
  
  # When offsets are not a vector of doubles
  offset <- rep("1", N)
  offset_2 <- NA
  offset_3 <- Inf
  offset_4 <- rep(NaN, N)
  
  # Then error thrown when trying to calculate the loss
  expect_error(loss(resps, preds, weights, offset, dist))
  expect_error(loss(resps, preds, weights, offset_2, dist))
  expect_error(loss(resps, preds, weights, offset_3, dist))
  expect_error(loss(resps, preds, weights, offset_4, dist))
})
test_that("Error thrown when baseline not vector of doubles", {
  # Given responses, weights, offset, predictions and default dist
  N <- 100
  dist <- gbm_dist()
  resps <- runif(N)
  weights <- runif(N)
  preds <- runif(N)
  offset <- runif(N)
  
  # When baseline are not a vector of doubles
  baseline <- rep("1", N)
  baseline_2 <- NA
  baseline_3 <- Inf
  baseline_4 <- rep(NaN, N)
  
  # Then error thrown when trying to calculate the loss
  expect_error(loss(resps, preds, weights, offset, dist, baseline))
  expect_error(loss(resps, preds, weights, offset, dist, baseline_2))
  expect_error(loss(resps, preds, weights, offset, dist, baseline_3))
  expect_error(loss(resps, preds, weights, offset, dist, baseline_4))
})
test_that("Error thrown when distribution is not a GBMDist", {
  # Given responses, weights, offset, predictions and default dist
  N <- 100
  dist <- gbm_dist()
  resps <- runif(N)
  weights <- runif(N)
  preds <- runif(N)
  offset <- runif(N)
  
  # When dist obj is not GBMDist
  class(dist) <- ""
  
  # Then error thrown when trying to calculate the loss
  expect_error(loss(resps, preds, weights, offset, dist))
})
test_that("Error thrown when length of responses, predictions, weights and baseline not the same", {
  # Given responses, weights, offset, predictions and a baselines
  # of varying lengths
  N <- 100
  N2 <- 102
  dist <- gbm_dist()
  resps <- runif(N)
  weights <- runif(N)
  preds <- runif(N)
  offset <- runif(N)
  baseline <- runif(N)
  
  resps_2 <- runif(N2)
  weights_2 <- runif(N2)
  preds_2 <- runif(N2)
  offset_2 <- runif(N2)
  baseline_2 <- runif(N2)
  
  # Then error thrown when trying to calculate the loss if 
  # length of responses, predictions, weights and baselines not the same
  expect_error(loss(resps_2, preds, weights, offset, dist, baseline))
  expect_error(loss(resps, preds_2, weights, offset, dist, baseline))
  expect_error(loss(resps, preds, weights_2, offset, dist, baseline))
  expect_error(loss(resps, preds_2, weights, offset_2, dist, baseline))
  expect_error(loss(resps, preds, weights, offset, dist, baseline_2))
})
test_that("Error thrown when length of offset is not the same as predictions", {
  # Given responses, weights, predictions and default dist
  N <- 100
  dist <- gbm_dist()
  resps <- runif(N)
  weights <- runif(N)
  preds <- runif(N)
  
  # When offset has different number of elements to preds
  offset <- runif(N+1)
  
  # Then error thrown when trying to calculate the loss
  expect_error(loss(resps, preds, weights, offset, dist))
})
test_that("Error thrown when distribution object not recognised - default method", {
  # Given responses, weights, predictions and offset
  N <- 100
  resps <- runif(N)
  weights <- runif(N)
  preds <- runif(N)
  offset <- runif(N)
  
  # When distribution not recognised
  dist <- gbm_dist()
  class(dist) <- c("WeirdGBMDist", "GBMDist")
  
  # Then error thrown when trying to calculate the loss
  expect_error(loss(resps, preds, weights, offset, dist))
})
test_that("Error thrown when distribution is CoxPH", {
  # Given responses, weights, predictions, offset and
  # CoxPH dist
  N <- 100
  resps <- runif(N)
  weights <- runif(N)
  preds <- runif(N)
  offset <- runif(N)
  dist <- gbm_dist("CoxPH")
  
  # Then error thrown when trying to calculate the loss
  expect_error(loss(resps, preds, weights, offset, dist))
})
test_that("Error thrown when distribution is Gamma", {
  # Given responses, weights, predictions, offset and
  # Gamma dist
  N <- 100
  resps <- runif(N)
  weights <- runif(N)
  preds <- runif(N)
  offset <- runif(N)
  dist <- gbm_dist("Gamma")
  
  # Then error thrown when trying to calculate the loss
  expect_error(loss(resps, preds, weights, offset, dist))
})
test_that("Error thrown when distribution is Huberized", {
  # Given responses, weights, predictions, offset and
  # Huberized dist
  N <- 100
  resps <- runif(N)
  weights <- runif(N)
  preds <- runif(N)
  offset <- runif(N)
  dist <- gbm_dist("Huberized")
  
  # Then error thrown when trying to calculate the loss
  expect_error(loss(resps, preds, weights, offset, dist))
})
test_that("Error thrown when distribution is Quantile", {
  # Given responses, weights, predictions, offset and
  # Quantile dist
  N <- 100
  resps <- runif(N)
  weights <- runif(N)
  preds <- runif(N)
  offset <- runif(N)
  dist <- gbm_dist("Quantile")
  
  # Then error thrown when trying to calculate the loss
  expect_error(loss(resps, preds, weights, offset, dist))
})
test_that("Error thrown when distribution is TDist", {
  # Given responses, weights, predictions, offset and
  # TDist dist
  N <- 100
  resps <- runif(N)
  weights <- runif(N)
  preds <- runif(N)
  offset <- runif(N)
  dist <- gbm_dist("TDist")
  
  # Then error thrown when trying to calculate the loss
  expect_error(loss(resps, preds, weights, offset, dist))
})
test_that("Error thrown when distribution is Tweedie", {
  # Given responses, weights, predictions, offset and
  # Tweedie dist
  N <- 100
  resps <- runif(N)
  weights <- runif(N)
  preds <- runif(N)
  offset <- runif(N)
  dist <- gbm_dist("Tweedie")
  
  # Then error thrown when trying to calculate the loss
  expect_error(loss(resps, preds, weights, offset, dist))
})

context("Check loss calculation correct for various distributions")
test_that("Correctly calculates AdaBoost loss", {
  # Given responses, weights, predictions, offset, baseline and
  # AdaBoost dist
  N <- 100
  resps <- runif(N)
  weights <- runif(N)
  preds <- runif(N)
  offset <- runif(N)
  baseline <- runif(N)
  dist <- gbm_dist("AdaBoost")
  
  # When calculting loss
  calc_loss <- loss(resps, preds, weights, offset, dist, baseline)
  
  # Then it is correct
  preds <- preds + offset
  loss_true <- weighted.mean(exp(-(2*resps-1)*preds), weights) - baseline
  expect_equal(calc_loss, loss_true)
})
test_that("Correctly calculates Bernoulli loss", {
  # Given responses, weights, predictions, offset, baseline and
  # AdaBoost dist
  N <- 100
  resps <- runif(N)
  weights <- runif(N)
  preds <- runif(N)
  offset <- runif(N)
  baseline <- runif(N)
  dist <- gbm_dist("AdaBoost")
  
  # When calculting loss
  calc_loss <- loss(resps, preds, weights, offset, dist, baseline)
  
  # Then it is correct
  preds <- preds + offset
  loss_true <- -2*weighted.mean(resps*preds - log(1+exp(preds)), weights) - baseline
  expect_equal(calc_loss, loss_true)
})

test_that("Correctly calculates Gaussian loss", {
  # Given responses, weights, predictions, offset, baseline and
  # AdaBoost dist
  N <- 100
  resps <- runif(N)
  weights <- runif(N)
  preds <- runif(N)
  offset <- runif(N)
  baseline <- runif(N)
  dist <- gbm_dist("AdaBoost")
  
  # When calculting loss
  calc_loss <- loss(resps, preds, weights, offset, dist, baseline)
  
  # Then it is correct
  preds <- preds + offset
  loss_true <- weighted.mean((resps - preds)^2, weights) - baseline
  expect_equal(calc_loss, loss_true)
})

test_that("Correctly calculates Laplace loss", {
  # Given responses, weights, predictions, offset, baseline and
  # AdaBoost dist
  N <- 100
  resps <- runif(N)
  weights <- runif(N)
  preds <- runif(N)
  offset <- runif(N)
  baseline <- runif(N)
  dist <- gbm_dist("AdaBoost")
  
  # When calculting loss
  calc_loss <- loss(resps, preds, weights, offset, dist, baseline)
  
  # Then it is correct
  preds <- preds + offset
  loss_true <- weighted.mean(abs(resps-preds), weights) - baseline
  expect_equal(calc_loss, loss_true)
})

test_that("Correctly calculates Poisson loss", {
  # Given responses, weights, predictions, offset, baseline and
  # AdaBoost dist
  N <- 100
  resps <- runif(N)
  weights <- runif(N)
  preds <- runif(N)
  offset <- runif(N)
  baseline <- runif(N)
  dist <- gbm_dist("AdaBoost")
  
  # When calculting loss
  calc_loss <- loss(resps, preds, weights, offset, dist, baseline)
  
  # Then it is correct
  preds <- preds + offset
  loss_true <- -2*weighted.mean(resps*preds-exp(preds), weights) - baseline
  expect_equal(calc_loss, loss_true)
})

