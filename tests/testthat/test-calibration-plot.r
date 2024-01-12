####################
# Author: James Hickey
#
# Series of tests to check if calibration plots working
#
####################


context("Testing calibration plot")
test_that("Error thrown if neither the knots or df specified - both NULL", {
  # Given data - based on example - but y and p not same length now
  kyph_dat <- data.frame(Kyphosis=sample(as.factor(c("absent", "present")), 100, replace=TRUE),
                         Age = sample(seq_len(175), 100, replace=TRUE), ncol=2)
  y <- as.numeric(kyph_dat$Kyphosis)-1
  x <- kyph_dat$Age
  glm1 <- glm(y~poly(x,2),family=binomial)
  p <- predict(glm1, type="response")
  
  # Then error thrown when both knots and df are NULL
  expect_error(calibrate_plot(y, p, df=NULL, knots=NULL, xlim=c(0,0.6), ylim=c(0,0.6)))
})
test_that("Error thrown if df is not a positive integer (if vector first element must be a positive integer)", {
  # Given data - based on example - but y and p not same length now
  kyph_dat <- data.frame(Kyphosis=sample(as.factor(c("absent", "present")), 100, replace=TRUE),
                         Age = sample(seq_len(175), 100, replace=TRUE), ncol=2)
  y <- as.numeric(kyph_dat$Kyphosis)-1
  x <- kyph_dat$Age
  glm1 <- glm(y~poly(x,2),family=binomial)
  p <- predict(glm1, type="response")
  
  # Then error thrown when df Not a positive integer
  expect_error(calibrate_plot(y, p, df=0, knots=NULL, xlim=c(0,0.6), ylim=c(0,0.6)))
  expect_error(calibrate_plot(y, p, df=1.4, knots=NULL, xlim=c(0,0.6), ylim=c(0,0.6)))
  expect_error(calibrate_plot(y, p, df="Wrong", knots=NULL, xlim=c(0,0.6), ylim=c(0,0.6)))
  expect_warning(calibrate_plot(y, p, df=c(1, 2), knots=NULL, xlim=c(0,0.6), ylim=c(0,0.6)))
  expect_error(calibrate_plot(y, p, df=c(1.4, 2), knots=NULL, xlim=c(0,0.6), ylim=c(0,0.6)))
  expect_error(calibrate_plot(y, p, df=NaN, knots=NULL, xlim=c(0,0.6), ylim=c(0,0.6)))
  expect_error(calibrate_plot(y, p, df=NA, knots=NULL, xlim=c(0,0.6), ylim=c(0,0.6)))
  expect_error(calibrate_plot(y, p, df=Inf, knots=NULL, xlim=c(0,0.6), ylim=c(0,0.6)))
})
test_that("Error thrown if y and p not same length", {
  # Given data - based on example - but y and p not same length now
  kyph_dat <- data.frame(Kyphosis=sample(as.factor(c("absent", "present")), 100, replace=TRUE),
                         Age = sample(seq_len(175), 100, replace=TRUE), ncol=2)
  y <- as.numeric(kyph_dat$Kyphosis)-1
  x <- kyph_dat$Age
  glm1 <- glm(y~poly(x,2),family=binomial)
  p <- predict(glm1, type="response")
  p <- p[seq(length(p)-1)]
  
  # Then error thrown
  expect_error(calibrate_plot(y, p))
})
test_that("Can run with defaults", {
  # Given data - based on example - but y and p not same length now
  kyph_dat <- data.frame(Kyphosis=sample(as.factor(c("absent", "present")), 100, replace=TRUE),
                         Age = sample(seq_len(175), 100, replace=TRUE), ncol=2)
  y <- as.numeric(kyph_dat$Kyphosis)-1
  x <- kyph_dat$Age
  glm1 <- glm(y~poly(x,2),family=binomial)
  p <- predict(glm1, type="response")
  
  
  # Then no error thrown
  expect_error(calibrate_plot(y, p), NA)
})
test_that("Can run with shade_col not NA", {
  # Given data - based on example - but y and p not same length now
  kyph_dat <- data.frame(Kyphosis=sample(as.factor(c("absent", "present")), 100, replace=TRUE),
                         Age = sample(seq_len(175), 100, replace=TRUE), ncol=2)
  y <- as.numeric(kyph_dat$Kyphosis)-1
  x <- kyph_dat$Age
  glm1 <- glm(y~poly(x,2),family=binomial)
  p <- predict(glm1, type="response")
  
  
  # Then no error thrown
  expect_error(calibrate_plot(y, p, shade_col=1), NA)
})
test_that("Can run with replace = FALSE", {
  # Given data - based on example - but y and p not same length now
  kyph_dat <- data.frame(Kyphosis=sample(as.factor(c("absent", "present")), 100, replace=TRUE),
                         Age = sample(seq_len(175), 100, replace=TRUE), ncol=2)
  y <- as.numeric(kyph_dat$Kyphosis)-1
  x <- kyph_dat$Age
  glm1 <- glm(y~poly(x,2),family=binomial)
  p <- predict(glm1, type="response")
  
  
  # Then no error thrown
  expect_error(calibrate_plot(y, p, replace=FALSE), NA)
})
test_that("Can run with shade_density != NULL", {
  # Given data - based on example - but y and p not same length now
  kyph_dat <- data.frame(Kyphosis=sample(as.factor(c("absent", "present")), 100, replace=TRUE),
                         Age = sample(seq_len(175), 100, replace=TRUE), ncol=2)
  y <- as.numeric(kyph_dat$Kyphosis)-1
  x <- kyph_dat$Age
  glm1 <- glm(y~poly(x,2),family=binomial)
  p <- predict(glm1, type="response")
  
  
  # Then no error thrown
  expect_error(calibrate_plot(y, p, shade_density=2.0), NA)
})
test_that("Can run  all distributions", {
  # Given data - based on example - but y and p not same length now
  kyph_dat <- data.frame(Kyphosis=sample(as.factor(c("absent", "present")), 100, replace=TRUE),
                         Age = sample(seq_len(175), 100, replace=TRUE), ncol=2)
  y <- as.numeric(kyph_dat$Kyphosis)-1
  x <- kyph_dat$Age
  glm1 <- glm(y~poly(x,2),family=binomial)
  p <- predict(glm1, type="response")
  
  
  # Then no error thrown
  expect_error(calibrate_plot(y, p, distribution = "AdaBoost"), NA)
  expect_error(calibrate_plot(y, p, distribution = "Bernoulli"), NA)
  expect_error(calibrate_plot(y, p, distribution = "CoxPH"), NA)
  expect_error(calibrate_plot(y, p, distribution = "Gamma"), NA)
  expect_error(calibrate_plot(y, p, distribution = "Gaussian"), NA)
  expect_error(calibrate_plot(y, p, distribution = "Laplace"), NA)
  expect_error(calibrate_plot(y, p, distribution = "Huberized"), NA)
  expect_error(calibrate_plot(y, p, distribution = "Pairwise"), NA)
  expect_error(calibrate_plot(y, p, distribution = "Poisson"), NA)
  expect_error(calibrate_plot(y, p, distribution = "Quantile"), NA)
  expect_error(calibrate_plot(y, p, distribution = "TDist"), NA)
  expect_error(calibrate_plot(y, p, distribution = "Tweedie"), NA)
})