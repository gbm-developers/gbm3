####################
# Author: James Hickey
#
# Series of tests to check if quantile rug plot is running correctly
#
####################

context("Testing quantile_rug plot")
test_that("quantile_rug runs without error", {
  probs <- rnorm(100)
  probs <- abs(probs)/sum(abs(probs))
  
  expect_error(quantile_rug(rnorm(100)), NA)
  expect_error(quantile_rug(rnorm(100), probs), NA)
})
test_that("quantile_rug throws error if probabilities outsider [0, 1]", {
  probs <- rnorm(100)
  probs[1] <- 5
  expect_error(quantile_rug(rnorm(100), probs))
})
test_that("output from quantile_rug is correct - no jittering", {
  probs <- (0:10)/10
  x <- rnorm(100)
  output <- quantile_rug(x, probs)
  expect_equal(output, rug(quantile(x[!is.na(x)], probs)))
})
test_that("quantile_rug jitters the inputs if quantiles < length(probabilities)", {
  set.seed(1)
  x <- rep(1:5, 2)
  probs <- (0:100)/100
  output <- quantile_rug(x, probs)
  
  set.seed(1)
  true_output <- quantile(x, probs)
  true_output <- jitter(true_output)
  
  expect_equal(output, rug(true_output))
})

