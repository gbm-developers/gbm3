####################
# Author: James Hickey
#
# Series of tests to check creation of GBMDist object
# for old gbm API.
#
####################
context("Test the creation of GBMDist objects for old API")
test_that("Error thrown if distribution not recognised", {
  # Given a made up distribution
  dist <- list(name="woooooo! Error!")
  
  # When create_dist_obj_for_gbmt_fit is called
  # Then an error is thrown
  expect_error(create_dist_obj_for_gbmt_fit(dist), 
               paste("The distribution",dist$name,"is not available in the gbm package."))
})
test_that("Can create default distributions", {
  # Given default distributions
  # When creating similar defaults from create_dist_obj_for_gbmt_fit
  # Then successfully create those objects
  expect_equal(gbm_dist("AdaBoost"), create_dist_obj_for_gbmt_fit(list(name="AdaBoost")))
  expect_equal(gbm_dist("Bernoulli"), create_dist_obj_for_gbmt_fit(list(name="Bernoulli")))
  expect_equal(gbm_dist("CoxPH"), create_dist_obj_for_gbmt_fit(list(name="CoxPH")))
  expect_equal(gbm_dist("Gamma"), create_dist_obj_for_gbmt_fit(list(name="Gamma")))
  expect_equal(gbm_dist("Gaussian"), create_dist_obj_for_gbmt_fit(list(name="Gaussian")))
  expect_equal(gbm_dist("Huberized"), create_dist_obj_for_gbmt_fit(list(name="Huberized")))
  expect_equal(gbm_dist("Laplace"), create_dist_obj_for_gbmt_fit(list(name="Laplace")))
  expect_equal(gbm_dist("Pairwise"), create_dist_obj_for_gbmt_fit(list(name="Pairwise")))
  expect_equal(gbm_dist("Poisson"), create_dist_obj_for_gbmt_fit(list(name="Poisson")))
  expect_equal(gbm_dist("Quantile"), create_dist_obj_for_gbmt_fit(list(name="Quantile")))
  expect_equal(gbm_dist("TDist"), create_dist_obj_for_gbmt_fit(list(name="TDist")))
  expect_equal(gbm_dist("Tweedie"), create_dist_obj_for_gbmt_fit(list(name="Tweedie")))
})
test_that("Can specify model dependent data - CoxPH", {
  # Given a distribution list with strata, tied.times.method
  # and prior.node.coeff.var
  dist <- list(name="CoxPH")
  strata <- c(1, 2, 1)
  tied.times.method <- "breslow"
  prior.node.coeff.var <- 150.4
  
  # When distribution object is created
  dist_obj <- create_dist_obj_for_gbmt_fit(dist, tied.times.method, 
                                           strata, prior.node.coeff.var)
  
  # Then it has the correct model dependent data
  expect_equal(dist_obj$original_strata_id, strata)
  expect_equal(dist_obj$ties, tied.times.method)
  expect_equal(dist_obj$prior_node_coeff_var, prior.node.coeff.var)
})
test_that("Can specify model dependent data - Pairwise", {
  # Given a distribution list with metric, max.rank and group
  dist <- list(name="Pairwise", metric="mrr", max.rank=1, group="query")
  
  # When distribution object is created
  dist_obj <- create_dist_obj_for_gbmt_fit(dist)
  
  # Then it has the correct model dependent data
  expect_equal(dist_obj$metric, dist$metric)
  expect_equal(dist_obj$group, dist$group)
  expect_equal(dist_obj$max_rank, dist$max.rank)
})
test_that("Can specify model dependent data - Quantile", {
  # Given a distribution list with alpha specified
  dist <- list(name="Quantile", alpha=0.5)
  
  # When distribution object is created
  dist_obj <- create_dist_obj_for_gbmt_fit(dist)
  
  # Then it has the correct alpha
  expect_equal(dist_obj$alpha, dist$alpha)
})
test_that("Can specify model dependent data - TDist", {
  # Given a distribution list with df specified
  dist <- list(name="TDist", df=10)
  
  # When distribution object is created
  dist_obj <- create_dist_obj_for_gbmt_fit(dist)
  
  # Then it has the correct df
  expect_equal(dist_obj$df, dist$df)
})
test_that("Can specify model dependent data - Tweedie", {
  # Given a distribution list with power specified
  dist <- list(name="Tweedie", power=100)
  
  # When distribution object is created
  dist_obj <- create_dist_obj_for_gbmt_fit(dist)
  
  # Then it has the correct power
  expect_equal(dist_obj$power, dist$power)
})
test_that("Case of distribution name characters is irrelevant", {
  # Creating distributions - then case of letters in name doesn't 
  # matter
  expect_equal(create_dist_obj_for_gbmt_fit(list(name="adaboost")), 
               create_dist_obj_for_gbmt_fit(list(name="AdaBoost")))
  expect_equal(create_dist_obj_for_gbmt_fit(list(name="bernoulli")),
               create_dist_obj_for_gbmt_fit(list(name="Bernoulli")))
  expect_equal(create_dist_obj_for_gbmt_fit(list(name="coxph")),
               create_dist_obj_for_gbmt_fit(list(name="CoxPH")))
  expect_equal(create_dist_obj_for_gbmt_fit(list(name="gamma")),
               create_dist_obj_for_gbmt_fit(list(name="Gamma")))
  expect_equal(create_dist_obj_for_gbmt_fit(list(name="gaussian")),
               create_dist_obj_for_gbmt_fit(list(name="Gaussian")))
  expect_equal(create_dist_obj_for_gbmt_fit(list(name="huberized")),
               create_dist_obj_for_gbmt_fit(list(name="Huberized")))
  expect_equal(create_dist_obj_for_gbmt_fit(list(name="laplace")),
               create_dist_obj_for_gbmt_fit(list(name="Laplace")))
  expect_equal(create_dist_obj_for_gbmt_fit(list(name="pairwise")),
               create_dist_obj_for_gbmt_fit(list(name="Pairwise")))
  expect_equal(create_dist_obj_for_gbmt_fit(list(name="poisson")),
               create_dist_obj_for_gbmt_fit(list(name="Poisson")))
  expect_equal(create_dist_obj_for_gbmt_fit(list(name="quantile")),
               create_dist_obj_for_gbmt_fit(list(name="Quantile")))
  expect_equal(create_dist_obj_for_gbmt_fit(list(name="tdist")), 
               create_dist_obj_for_gbmt_fit(list(name="TDist")))
  expect_equal(create_dist_obj_for_gbmt_fit(list(name="tweedie")),
               create_dist_obj_for_gbmt_fit(list(name="Tweedie")))
})
