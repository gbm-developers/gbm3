context("Test available distributions")
test_that("available_distributions function returns vector of names of implemented distributions in package", {
  implemented_distributions <- c("AdaBoost", "Bernoulli", "CoxPH", "Gamma", "Gaussian", "Huberized", "Laplace", "Pairwise",
                                 "Poisson", "Quantile", "TDist", "Tweedie")
  expect_equal(available_distributions(), implemented_distributions)
})