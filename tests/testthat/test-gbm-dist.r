####################
# Author: James Hickey
#
# Series of test to validate the GBMDist objects.
#
####################

context("Testing GBMDist Object Definition:")

test_that("Default distribution is Gaussian", {
  gbm_dist_obj <- gbm_dist()
  expect_equal(gbm_dist_obj$name, "Gaussian")
})

test_that("Check AdaBoost Distribution Object has correct class attributes", {
  gbm_dist_obj <- gbm_dist(name="AdaBoost")
  expect_true("GBMDist" %in% class(gbm_dist_obj))
  expect_true("AdaBoostGBMDist" %in% class(gbm_dist_obj))
})

test_that("Check Bernoulli Distribution Object has correct class attributes", {
  gbm_dist_obj <- gbm_dist(name="Bernoulli")
  expect_true("GBMDist" %in% class(gbm_dist_obj))
  expect_true("BernoulliGBMDist" %in% class(gbm_dist_obj))
})

test_that("Check CoxPH Distribution Object has correct class attributes", {
  gbm_dist_obj <- gbm_dist(name="CoxPH")
  expect_true("GBMDist" %in% class(gbm_dist_obj))
  expect_true("CoxPHGBMDist" %in% class(gbm_dist_obj))
})

test_that("Check Gamma Distribution Object has correct class attributes", {
  gbm_dist_obj <- gbm_dist(name="Gamma")
  expect_true("GBMDist" %in% class(gbm_dist_obj))
  expect_true("GammaGBMDist" %in% class(gbm_dist_obj))
})

test_that("Check Gaussian Distribution Object has correct class attributes", {
  gbm_dist_obj <- gbm_dist(name="Gaussian")
  expect_true("GBMDist" %in% class(gbm_dist_obj))
  expect_true("GaussianGBMDist" %in% class(gbm_dist_obj))
})

test_that("Check Huberized Distribution Object has correct class attributes", {
  gbm_dist_obj <- gbm_dist(name="Huberized")
  expect_true("GBMDist" %in% class(gbm_dist_obj))
  expect_true("HuberizedGBMDist" %in% class(gbm_dist_obj))
})

test_that("Check Laplace Distribution Object has correct class attributes", {
  gbm_dist_obj <- gbm_dist(name="Laplace")
  expect_true("GBMDist" %in% class(gbm_dist_obj))
  expect_true("LaplaceGBMDist" %in% class(gbm_dist_obj))
})
test_that("Check Pairwise Distribution Object has correct class attributes", {
  gbm_dist_obj <- gbm_dist(name="Pairwise")
  expect_true("GBMDist" %in% class(gbm_dist_obj))
  expect_true("PairwiseGBMDist" %in% class(gbm_dist_obj))
})
test_that("Check Poisson Distribution Object has correct class attributes", {
  gbm_dist_obj <- gbm_dist(name="Poisson")
  expect_true("GBMDist" %in% class(gbm_dist_obj))
  expect_true("PoissonGBMDist" %in% class(gbm_dist_obj))
})
test_that("Check Quantile Distribution Object has correct class attributes", {
  gbm_dist_obj <- gbm_dist(name="Quantile")
  expect_true("GBMDist" %in% class(gbm_dist_obj))
  expect_true("QuantileGBMDist" %in% class(gbm_dist_obj))
})

test_that("Check TDist Distribution Object has correct class attributes", {
  gbm_dist_obj <- gbm_dist(name="TDist")
  expect_true("GBMDist" %in% class(gbm_dist_obj))
  expect_true("TDistGBMDist" %in% class(gbm_dist_obj))
})

test_that("Check Tweedie Distribution Object has correct class attributes", {
  gbm_dist_obj <- gbm_dist(name="Tweedie")
  expect_true("GBMDist" %in% class(gbm_dist_obj))
  expect_true("TweedieGBMDist" %in% class(gbm_dist_obj))
})