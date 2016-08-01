####################
# Author: James Hickey
#
# Series of tests to check if how training params are updated if groups are present
#
####################

context("Input checking update_num_train_groups")
test_that("update_num_train_groups throws an error if not given a GBMTrainParams object", {
  # Given a default train params and dist object
  params <- training_params()
  dist <- gbm_dist()
  
  # When training params stripped of class and update_num_train_groups called
  class(params) <- "wrong"
  
  # Then error thrown
  expect_error(update_num_train_groups(params, dist))
})
test_that("update_num_train_groups throws an error if not given a GBMDist object", {
  # Given a default train params and dist object
  params <- training_params()
  dist <- gbm_dist()
  
  # When dist is stripped of class and update_num_train_groups called
  class(params) <- "wrong"
  
  # Then error thrown
  expect_error(update_num_train_groups(params, dist))
})

context("Testing update_num_train_groups behaviour")
test_that("update_num_train_groups returns original params object when given an AdaBoostGBMDist obj", {
  # Given some training parameters and a default dist
  params <- training_params()
  dist <- gbm_dist("AdaBoost")
  
  # When training parameters are updated using update_num_train_groups
  # Then no update happens
  expect_equal(update_num_train_groups(params, dist), params)
})
test_that("update_num_train_groups returns original params object when given an BernoulliGBMDist obj", {
  # Given some training parameters and a default dist
  params <- training_params()
  dist <- gbm_dist("Bernoulli")
  
  # When training parameters are updated using update_num_train_groups
  # Then no update happens
  expect_equal(update_num_train_groups(params, dist), params)
})
test_that("update_num_train_groups returns original params object when given an CoxPHGBMDist obj", {
  # Given some training parameters and a default dist
  params <- training_params()
  dist <- gbm_dist("CoxPH")
  
  # When training parameters are updated using update_num_train_groups
  # Then no update happens
  expect_equal(update_num_train_groups(params, dist), params)
})
test_that("update_num_train_groups returns original params object when given an GammaGBMDist obj", {
  # Given some training parameters and a default dist
  params <- training_params()
  dist <- gbm_dist("Gamma")
  
  # When training parameters are updated using update_num_train_groups
  # Then no update happens
  expect_equal(update_num_train_groups(params, dist), params)
})
test_that("update_num_train_groups returns original params object when given an GaussianGBMDist obj", {
  # Given some training parameters and a default dist
  params <- training_params()
  dist <- gbm_dist("Gaussian")
  
  # When training parameters are updated using update_num_train_groups
  # Then no update happens
  expect_equal(update_num_train_groups(params, dist), params)
})
test_that("update_num_train_groups returns original params object when given an HuberizedGBMDist obj", {
  # Given some training parameters and a default dist
  params <- training_params()
  dist <- gbm_dist("Huberized")
  
  # When training parameters are updated using update_num_train_groups
  # Then no update happens
  expect_equal(update_num_train_groups(params, dist), params)
})
test_that("update_num_train_groups returns original params object when given an LaplaceGBMDist obj", {
  # Given some training parameters and a default dist
  params <- training_params()
  dist <- gbm_dist("Laplace")
  
  # When training parameters are updated using update_num_train_groups
  # Then no update happens
  expect_equal(update_num_train_groups(params, dist), params)
})
test_that("update_num_train_groups returns original params object when given an PoissonGBMDist obj", {
  # Given some training parameters and a default dist
  params <- training_params()
  dist <- gbm_dist("Poisson")
  
  # When training parameters are updated using update_num_train_groups
  # Then no update happens
  expect_equal(update_num_train_groups(params, dist), params)
})
test_that("update_num_train_groups returns original params object when given an QuantileGBMDist obj", {
  # Given some training parameters and a default dist
  params <- training_params()
  dist <- gbm_dist("Quantile")
  
  # When training parameters are updated using update_num_train_groups
  # Then no update happens
  expect_equal(update_num_train_groups(params, dist), params)
})
test_that("update_num_train_groups returns original params object when given an TDistGBMDist obj", {
  # Given some training parameters and a default dist
  params <- training_params()
  dist <- gbm_dist("TDist")
  
  # When training parameters are updated using update_num_train_groups
  # Then no update happens
  expect_equal(update_num_train_groups(params, dist), params)
})
test_that("update_num_train_groups returns original params object when given an TweedieGBMDist obj", {
  # Given some training parameters and a default dist
  params <- training_params()
  dist <- gbm_dist("AdaBoost")
  
  # When training parameters are updated using update_num_train_groups
  # Then no update happens
  expect_equal(update_num_train_groups(params, dist), params)
})
test_that("update_num_train_groups throws error if training_params id indicates multiple rows per observation - Pairwise", {
  # Given a training params object and an appropriately set-up pairwise object
  # such that observations have multiple rows
  N <- 100
  dist <- gbm_dist("Pairwise")
  params <- training_params(num_train=N, id=sample(seq_len(N/5), N, replace=TRUE))
  
  # When attempt to update_num_train_groups
  # Then error is thrown
  expect_error(update_num_train_groups(params, dist))
})
test_that("update_num_train_groups returns correctly updated training_params - Pairwise", {
  # Given some training parameters object and a Pairwise object with groups
  N <- 100
  params <- training_params(num_train=N/2, id=seq_len(N))
  dist <- gbm_dist("Pairwise")
  dist$groups <- sample(seq_len(5), N, replace=TRUE)
  
  # When update_num_train_groups
  params_update <- update_num_train_groups(params, dist)
  
  # Then correctly calculates update
  num_train_groups <- max(1, round(params$train_fraction * nlevels(dist$group)))
  params$num_train_rows <- max(which(dist$group == levels(dist$group)[num_train_groups]))
  params$num_train <- params$num_train_rows
  expect_equal(params_update, params)
})