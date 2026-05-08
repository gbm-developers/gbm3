####################
# Author: James Hickey
#
# Series of tests to test baseline hazard function
#
####################

context("Testing error checking - baseline_hazard")
test_that("error thrown if cumulative is not logical", {
  # Given inputs  - MADE UP
  N<-1
  surv_time <- runif(N)
  delta <- sample(0:1, N, replace=TRUE)
  cox_preds <- abs(runif(N))
  cox_preds <- cox_preds/(max(cox_preds) - min(cox_preds))
  
  # When cumulative is not logical 
  # Then error is thrown
  expect_error(baseline_hazard(surv_time, delta, cox_preds, cumulative = NA))
  expect_error(baseline_hazard(surv_time, delta, cox_preds, cumulative = NaN))
  expect_error(baseline_hazard(surv_time, delta, cox_preds, cumulative = "Wrong"))
  expect_error(baseline_hazard(surv_time, delta, cox_preds, cumulative = Inf))
  expect_error(baseline_hazard(surv_time, delta, cox_preds, cumulative = 1))
  expect_error(baseline_hazard(surv_time, delta, cox_preds, cumulative = c(TRUE, FALSE)))
})
test_that("error thrown if smooth is not logical", {
  # Given inputs  - MADE UP
  N<-1
  surv_time <- runif(N)
  delta <- sample(0:1, N, replace=TRUE)
  cox_preds <- abs(runif(N))
  cox_preds <- cox_preds/(max(cox_preds) - min(cox_preds))
  
  # When smooth is not logical 
  # Then error is thrown
  expect_error(baseline_hazard(surv_time, delta, cox_preds, smooth = NA))
  expect_error(baseline_hazard(surv_time, delta, cox_preds, smooth = NaN))
  expect_error(baseline_hazard(surv_time, delta, cox_preds, smooth = "Wrong"))
  expect_error(baseline_hazard(surv_time, delta, cox_preds, smooth = Inf))
  expect_error(baseline_hazard(surv_time, delta, cox_preds, smooth = 1))
  expect_error(baseline_hazard(surv_time, delta, cox_preds, smooth = c(TRUE, FALSE)))
})
test_that("error thrown if delta not in {0, 1}", {
  # Given inputs  - MADE UP
  N<-1
  surv_time <- runif(N)
  delta <- sample(0:1, N, replace=TRUE)
  cox_preds <- abs(runif(N))
  cox_preds <- cox_preds/(max(cox_preds) - min(cox_preds))
  
  # When all delta not in {0, 1}
  delta[1] <- 1.1
  # Then error is thrown
  expect_error(baseline_hazard(surv_time, delta, cox_preds))
})
test_that("error thrown if delta not an atomic", {
  # Given inputs  - MADE UP
  N<-1
  surv_time <- runif(N)
  delta <- sample(0:1, N, replace=TRUE)
  cox_preds <- abs(runif(N))
  cox_preds <- cox_preds/(max(cox_preds) - min(cox_preds))
  
  # When all delta not an atomic
  # Then error is thrown
  expect_error(baseline_hazard(surv_time, as.list(delta), cox_preds))
})
test_that("error thrown when smooth=FALSE, cumulative=FALSE and eval_times are not NULL", {
  # Given inputs  - MADE UP
  N<-1
  surv_time <- runif(N)
  delta <- sample(0:1, N, replace=TRUE)
  cox_preds <- abs(runif(N))
  cox_preds <- cox_preds/(max(cox_preds) - min(cox_preds))
  
  # When smooth=FALSE, cumulative=FALSE and eval_times not NULL 
  # Then error is thrown
  expect_error(baseline_hazard(surv_time, delta, cox_preds, eval_times = runif(N), 
                               cumulative = FALSE, smooth=FALSE))
})
test_that("error thrown when survival times are not atomic of doubles", {
  # Given inputs  - MADE UP
  N<-1
  surv_time <- runif(N)
  delta <- sample(0:1, N, replace=TRUE)
  cox_preds <- abs(runif(N))
  cox_preds <- cox_preds/(max(cox_preds) - min(cox_preds))
  
  # When survival_times not an atomic of doubls
  # Then error is thrown
  expect_error(baseline_hazard(as.list(surv_time), delta, cox_preds, eval_times = runif(N)))
  expect_error(baseline_hazard(Inf, delta, cox_preds, eval_times = runif(N)))
})

context("Testing it can run - baseline hazard")
test_that("can run cumulative=TRUE", {
  # Given inputs  - MADE UP
  N<-100
  surv_time <- runif(N)
  delta <- sample(0:1, N, replace=TRUE)
  cox_preds <- abs(runif(N))
  cox_preds <- cox_preds/(max(cox_preds) - min(cox_preds))
  
  # Whenrun cumulative=TRUE
  # Then no error is thrown
  expect_error(baseline_hazard(surv_time, delta, cox_preds, cumulative=TRUE), NA)
})

test_that("can run smooth=TRUE", {
  # Given inputs  - MADE UP
  N<-100
  surv_time <- runif(N)
  delta <- sample(0:1, N, replace=TRUE)
  cox_preds <- abs(runif(N))
  cox_preds <- cox_preds/(max(cox_preds) - min(cox_preds))
  
  # Whenrun cumulative=TRUE
  # Then no error is thrown
  expect_error(baseline_hazard(surv_time, delta, cox_preds, smooth=TRUE), NA)
})

test_that("can run with eval_times", {
  # Given inputs  - MADE UP
  N<-100
  surv_time <- runif(N)
  delta <- sample(0:1, N, replace=TRUE)
  cox_preds <- abs(runif(N))
  cox_preds <- cox_preds/(max(cox_preds) - min(cox_preds))
  
  # Whenrun cumulative=TRUE
  # Then no error is thrown
  expect_error(baseline_hazard(surv_time, delta, cox_preds, eval_times=surv_time+1), NA)
})

test_that("eval_times controls output locations", {
  surv_time <- c(1, 2, 3, 4)
  delta <- c(1, 1, 0, 1)
  cox_preds <- rep(0, length(surv_time))
  eval_times <- c(1.5, 3.5)

  hazard <- baseline_hazard(surv_time, delta, cox_preds,
                            eval_times=eval_times)

  expect_length(hazard, length(eval_times))
  expect_equal(hazard, c(5 / 12, 4 / 3))
})
