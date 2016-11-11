####################
# Author: James Hickey
#
# Series of tests to check if group determination 
# and weight updates based on groups
#
####################

context("Testing determine groups")
test_that("determine_groups throws error if not passed a GBMDist object", {
  # Given data and distribution obj
  # create query groups, with an average size of 25 items each
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
  
  data <- data.frame(Y, query=query, X1, X2, X3)
  dist <- gbm_dist("Pairwise", metric="ndcg", group="query")
  
  # When dist is stripped of class
  class(dist) <- "Wrong"
  
  # Then determine_groups throws an error
  expect_error(determine_groups(data, Y, dist))
})
test_that("determine_groups returns original distribution_obj if NOT Pairwise", {
  # Given data and distribution objs - NOT Pairwise
  # create query groups, with an average size of 25 items each
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
  
  data <- data.frame(Y, query=query, X1, X2, X3)
  dist_1 <- gbm_dist("AdaBoost")
  dist_2 <- gbm_dist("Bernoulli")
  dist_3 <- gbm_dist("CoxPH")
  dist_4 <- gbm_dist("Gamma")
  dist_5 <- gbm_dist("Gaussian")
  dist_6 <- gbm_dist("Huberized")
  dist_7 <- gbm_dist("Laplace")
  dist_8 <- gbm_dist("Poisson")
  dist_9 <- gbm_dist("Quantile")
  dist_10 <- gbm_dist("TDist")
  dist_11 <- gbm_dist("Tweedie")

  # When determine_groups called
  # Then returns original distribution
  expect_equal(determine_groups(data, Y, dist_1), dist_1)
  expect_equal(determine_groups(data, Y, dist_2), dist_2)
  expect_equal(determine_groups(data, Y, dist_3), dist_3)
  expect_equal(determine_groups(data, Y, dist_4), dist_4)
  expect_equal(determine_groups(data, Y, dist_5), dist_5)
  expect_equal(determine_groups(data, Y, dist_6), dist_6)
  expect_equal(determine_groups(data, Y, dist_7), dist_7)
  expect_equal(determine_groups(data, Y, dist_8), dist_8)
  expect_equal(determine_groups(data, Y, dist_9), dist_9)
  expect_equal(determine_groups(data, Y, dist_10), dist_10)
  expect_equal(determine_groups(data, Y, dist_11), dist_11)
})
test_that("determine_groups throws error if group field in distribution_obj is NULL - Pairwise", {
  # Given data and distribution obj - Pairwise (group set to NULL)
  # create query groups, with an average size of 25 items each
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
  
  data <- data.frame(Y, query=query, X1, X2, X3)
  dist <- gbm_dist("Pairwise", metric="ndcg", group="query")
  dist$group <- NULL
  
  # When determine_groups called
  # Then determine_groups throws an error
  expect_error(determine_groups(data, Y, dist))
})
test_that("determine_groups throws an error if group does not occur as column in original_data - Pairwise", {
  # Given data and distribution obj - Pairwise with group not in original_data
  # create query groups, with an average size of 25 items each
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
  
  data <- data.frame(Y, query=query, X1, X2, X3)
  dist <- gbm_dist("Pairwise", metric="ndcg", group="WRONG")
  
  # When determine_groups is called
  # Then determine_groups throws an error
  expect_error(determine_groups(data, Y, dist))
})
test_that("determine_groups can run and populates group_index & group_order fields and updates group field - Pairwise", {
  # Given data and distribution obj - Pairwise and group correct
  # create query groups, with an average size of 25 items each
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
  
  data <- data.frame(Y, query=query, X1, X2, X3)
  dist <- gbm_dist("Pairwise", metric="ndcg", group="query")
  
  # When determine_groups called
  updated_dist <- determine_groups(data, Y, dist)
  
  # Then group_index and group_order fields created and group updated
  expect_false(is.null(updated_dist$group_index))
  expect_false(is.null(updated_dist$group_order))
  expect_false(any(updated_dist$group == dist$group))
})

context("Testing weight_group_consistency")
test_that("weight_group_consistency throws an error if distribution object is not GBMDist class", {
  # Given weights and a distribution_object
  w <- rep(1, 10)
  dist <- gbm_dist()
  
  # When dist obj has class removed
  class(dist) <- "Wrong"
  
  # Then weight_group_consistency throws error
  expect_error(weight_group_consistency(dist, w))
})
test_that("weight_group_consistency does not alter weights if NOT Pairwise", {
  # Given weights and distribution objects NOT Pairwise
  w <- rep(1, 10)
  dist_1 <- gbm_dist("AdaBoost")
  dist_2 <- gbm_dist("Bernoulli")
  dist_3 <- gbm_dist("CoxPH")
  dist_4 <- gbm_dist("Gamma")
  dist_5 <- gbm_dist("Gaussian")
  dist_6 <- gbm_dist("Huberized")
  dist_7 <- gbm_dist("Laplace")
  dist_8 <- gbm_dist("Poisson")
  dist_9 <- gbm_dist("Quantile")
  dist_10 <- gbm_dist("TDist")
  dist_11 <- gbm_dist("Tweedie")
  
  # When weight_group_consistency is called
  # Then the original weights are returned
  expect_equal(weight_group_consistency(dist_1, w), w)
  expect_equal(weight_group_consistency(dist_2, w), w)
  expect_equal(weight_group_consistency(dist_3, w), w)
  expect_equal(weight_group_consistency(dist_4, w), w)
  expect_equal(weight_group_consistency(dist_5, w), w)
  expect_equal(weight_group_consistency(dist_6, w), w)
  expect_equal(weight_group_consistency(dist_7, w), w)
  expect_equal(weight_group_consistency(dist_8, w), w)
  expect_equal(weight_group_consistency(dist_9, w), w)
  expect_equal(weight_group_consistency(dist_10, w), w)
  expect_equal(weight_group_consistency(dist_11, w), w)
})
test_that("weight_group_consistency throws error if weights across groups are different - Pairwise", {
  # Given weights are different across groups and a Pairwise distribution
  # create query groups, with an average size of 25 items each
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
  data <- data.frame(Y, query=query, X1, X2, X3)
  dist <- gbm_dist("Pairwise", metric="ndcg", group="query")
  dist <- determine_groups(data, Y, dist)
  
  # When weights are not same across groups
  w <- rep(1, N)
  w[dist$group == 1] <- abs(rnorm(length(dist$group[dist$group == 1])))
  
  # Then error is thrown
  expect_error(weight_group_consistency(dist, w))
})

