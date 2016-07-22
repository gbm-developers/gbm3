####################
# Author: James Hickey
#
# Series of tests for get_misc functionality
#
####################

context("Testing ability to get misc for a distribution object")
test_that("Error thrown when passed object not GBMDist class", {
  # Given a distribution object
  dist <- gbm_dist()
  
  # When it has its class removed
  class(dist) <- ""
  
  # Then error is thrown on get_misc
  expect_error(get_misc(dist))
})
test_that("Error thrown for unrecognised distribution object", {
  # Given a distribution object
  dist <- gbm_dist()
  
  # When it has an unrecognised class
  class(dist) <- c("WeirdGBMDist", "GBMDist")
  
  # Then error is thrown on get_misc
  expect_error(get_misc(dist))
})
test_that("get_misc returns a list", {
  # Given a distribution object
  # for each distribution
  dist_1 <- gbm_dist("AdaBoost")
  dist_2 <- gbm_dist("Bernoulli")
  dist_3 <- gbm_dist("CoxPH")
  dist_4 <- gbm_dist("Gamma")
  dist_5 <- gbm_dist("Gaussian")
  dist_6 <- gbm_dist("Huberized")
  dist_7 <- gbm_dist("Laplace")
  dist_8 <- gbm_dist("Pairwise")
  dist_9 <- gbm_dist("Poisson")
  dist_10 <- gbm_dist("Quantile")
  dist_11 <- gbm_dist("TDist")
  dist_12 <- gbm_dist("Tweedie")
    
  # When it misc is gotten
  # Then returns a list
  expect_true(is.list(get_misc(dist_1)))
  expect_true(is.list(get_misc(dist_2)))
  expect_true(is.list(get_misc(dist_3)))
  expect_true(is.list(get_misc(dist_4)))
  expect_true(is.list(get_misc(dist_5)))
  expect_true(is.list(get_misc(dist_6)))
  expect_true(is.list(get_misc(dist_7)))
  expect_true(is.list(get_misc(dist_8)))
  expect_true(is.list(get_misc(dist_9)))
  expect_true(is.list(get_misc(dist_10)))
  expect_true(is.list(get_misc(dist_11)))
  expect_true(is.list(get_misc(dist_12)))
  
})
test_that("Can get misc for AdaBoost", {
  # Given an AdaBoost dist
  dist <- gbm_dist("AdaBoost")
  
  # Then get_misc finds NA
  expect_true(is.na(get_misc(dist)))
})
test_that("Can get misc for Bernoulli", {
  # Given an Bernoulli dist
  dist <- gbm_dist("Bernoulli")
  
  # Then get_misc finds NA
  expect_true(is.na(get_misc(dist)))
})
test_that("Can get misc for CoxPH", {
  # Given a CoXPH distribution
  dist <- gbm_dist("CoxPH", ties="breslow")
  
  # When get misc
  misc <- get_misc(dist)
  
  # Then it is the ties
  expect_equal(misc$ties, "breslow")
})
test_that("Can get misc for Gamma", {
  # Given an Gamma dist
  dist <- gbm_dist("Gamma")
  
  # Then get_misc finds NA
  expect_true(is.na(get_misc(dist)))
})
test_that("Can get misc for Gaussian", {
  # Given an Gaussian dist
  dist <- gbm_dist("Gaussian")
  
  # Then get_misc finds NA
  expect_true(is.na(get_misc(dist)))
})
test_that("Can get misc for Huberized", {
  # Given an Huberized dist
  dist <- gbm_dist("Huberized")
  
  # Then get_misc finds NA
  expect_true(is.na(get_misc(dist)))
})
test_that("Can get misc for Laplace", {
  # Given an Laplace dist
  dist <- gbm_dist("Laplace")
  
  # Then get_misc finds NA
  expect_true(is.na(get_misc(dist)))
})
test_that("Can get misc for Pairwise", {
  # Given a pairwise dist
  group <- "query"
  max_rank <- 2
  dist <- gbm_dist("Pairwise", group="query", max.rank=max_rank)
  
  # When misc is gotten
  misc <- get_misc(dist)
  
  # Then it is a c(group, max_rank) in a list
  expect_equal(misc$GroupsAndRanks, c(group, max_rank))
})
test_that("Can get misc for Poisson", {
  # Given an Poisson dist
  dist <- gbm_dist("Poisson")
  
  # Then get_misc finds NA
  expect_true(is.na(get_misc(dist)))
})
test_that("Can get misc for Quantile", {
  # Given a Quantile dist
  dist <- gbm_dist("Quantile", alpha=0.25)
  
  # When misc gotten
  misc <- get_misc(dist)
  
  # Then misc is alpha
  expect_equal(misc$alpha, 0.25)
  
})
test_that("Can get misc for TDist", {
  # Given a TDist dist
  df <- 4
  dist <- gbm_dist("TDist", df=df)
  
  # When misc is gotten
  misc <- get_misc(dist)
  
  # Then it is df
  expect_equal(misc$df, df)
})
test_that("Can get misc for Tweedie", {
  # Given a Tweedie distribution
  power <- 1.7
  dist <- gbm_dist("Tweedie", power=power)
  
  # When misc is gotten
  misc <- get_misc(dist)
  
  # Then it is the power of the dist
  expect_equal(misc$power, power)
})
