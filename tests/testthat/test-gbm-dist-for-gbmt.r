####################
# Author: James Hickey
#
# Series of tests to check the creation of distribution objects in gbm
#
####################

context("Testing creation of dist objects in gbm")
test_that("Error thrown if name in distribution is null", {
  dist <- list()
  expect_error(create_dist_obj_for_gbmt_fit(dist))
})
test_that("Error thrown if name not recognised", {
  dist <- list(name="Wrong")
  expect_error(create_dist_obj_for_gbmt_fit(dist))
})
test_that("Can create default CoxPH object", {
  def_cox_ph <- gbm_dist("CoxPH")
  dist <- create_dist_obj_for_gbmt_fit(list(name="CoxPH"))
  expect_equal(def_cox_ph, dist)
})
test_that("Can create default TDist object", {
  def_tdist_ph <- gbm_dist("TDist")
  dist <- create_dist_obj_for_gbmt_fit(list(name="TDist"))
  expect_equal(def_tdist_ph, dist)
})
test_that("Can create default Pairwise object", {
  def_pairwise_ph <- gbm_dist("Pairwise")
  dist <- create_dist_obj_for_gbmt_fit(list(name="Pairwise"))
  expect_equal(def_pairwise_ph, dist)
})
test_that("Can create default Quantile object", {
  def_quantile_ph <- gbm_dist("Quantile")
  dist <- create_dist_obj_for_gbmt_fit(list(name="Quantile"))
  expect_equal(def_quantile_ph, dist)
})
test_that("Can create default Tweedie object", {
  def_tweedie_ph <- gbm_dist("Tweedie")
  dist <- create_dist_obj_for_gbmt_fit(list(name="Tweedie"))
  expect_equal(def_tweedie_ph, dist)
})
test_that("Can create other distributions", {
  expect_equal(create_dist_obj_for_gbmt_fit(list(name="AdaBoost")), gbm_dist("AdaBoost"))
  expect_equal(create_dist_obj_for_gbmt_fit(list(name="Bernoulli")), gbm_dist("Bernoulli"))
  expect_equal(create_dist_obj_for_gbmt_fit(list(name="Gamma")), gbm_dist("Gamma"))
  expect_equal(create_dist_obj_for_gbmt_fit(list(name="Gaussian")), gbm_dist("Gaussian"))
  expect_equal(create_dist_obj_for_gbmt_fit(list(name="Laplace")), gbm_dist("Laplace"))
  expect_equal(create_dist_obj_for_gbmt_fit(list(name="Huberized")), gbm_dist("Huberized"))
  expect_equal(create_dist_obj_for_gbmt_fit(list(name="Poisson")), gbm_dist("Poisson"))
})
test_that("Can set misc", {
  expect_equal(create_dist_obj_for_gbmt_fit(list(name="Quantile", alpha=1.0)), gbm_dist("Quantile", alpha=1.0))
  expect_equal(create_dist_obj_for_gbmt_fit(list(name="TDist", df=4.5)), gbm_dist("TDist", df=4.5))
  expect_equal(create_dist_obj_for_gbmt_fit(list(name="Tweedie", power=3.0)), gbm_dist("Tweedie", power=3.0))
  expect_equal(create_dist_obj_for_gbmt_fit(list(name="Pairwise", metric="ndcg", max.rank=1, group="g")),
               gbm_dist("Pairwise", metric="ndcg", max_rank=1, group="g"))
  expect_equal(create_dist_obj_for_gbmt_fit(list(name="CoxPH"), tied.times.method="breslow", prior.node.coeff.var=10, strata=c(1, 2)),
               gbm_dist("CoxPH", ties="breslow", strata=c(1,2), prior_node_coeff_var=10))
})



  