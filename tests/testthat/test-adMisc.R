##############
# Author: James Hickey
# 
# Series of tests that flesh out the behaviour of the parameter misc.
# Help identify if the refactoring is changing the higher level behaviour.
#
##############
context("Testing adMisc:")
test_that("The misc parameter does not impact gbm fitting with - Adaboost distribution", {
  # create some data
  set.seed(1)
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  
  # random weights if you want to experiment with them
  w <- rexp(N)
  w <- N*w/sum(w)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  
  # Fit gbm with and without misc
  admisc <- runif(N)
  set.seed(15) 
  gbm.no.misc <- gbm.fit(x=data[,-1], y=data[,1], distribution = "adaboost")
  
  set.seed(15)
  gbm.with.misc <- gbm.fit(x=data[,-1], y=data[,1], misc=admisc, distribution = "adaboost")
  expect_equal(gbm.no.misc, gbm.with.misc)
})
test_that("The misc parameter does not impact gbm fitting with - Bernoulli distribution", {
  # create some data
  set.seed(1)
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  
  # random weights if you want to experiment with them
  w <- rexp(N)
  w <- N*w/sum(w)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  
  # Fit gbm with and without misc
  admisc <- runif(N)
  set.seed(15) 
  gbm.no.misc <- gbm.fit(x=data[,-1], y=data[,1], distribution = "bernoulli")
  
  set.seed(15)
  gbm.with.misc <- gbm.fit(x=data[,-1], y=data[,1], misc=admisc, distribution = "bernoulli")
  expect_equal(gbm.no.misc, gbm.with.misc)
})
test_that("The misc parameter does not impact gbm fitting when it passes censoring indicator - Cox Partial Hazards model", {
  require(survival)
  
  ## GIVEN
  # create some data
  N <- 3000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  f <- 0.5*sin(3*X1 + 5*X2^2 + mu/10)
  tt.surv <- rexp(N,exp(f))
  tt.cens <- rexp(N,0.5)
  delta <- as.numeric(tt.surv <= tt.cens)
  tt <- apply(cbind(tt.surv,tt.cens),1,min)
  
  # throw in some missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  # random weights if you want to experiment with them
  w <- rep(1,N)
  
  data <- data.frame(y=Surv(tt,delta),X1=X1,X2=X2,X3=X3)
  
  ## WHEN
  # Fit gbm with and without misc
  set.seed(15) 
  gbm.no.misc <- gbm.fit(x=data[,-1], y=data[,1], distribution = "coxph")
  
  set.seed(15)
  gbm.with.misc <- gbm.fit(x=data[,-1], y=data[,1], misc=delta, distribution = "coxph")
  
  ## THEN
  expect_true(isTRUE(all.equal(gbm.no.misc, gbm.with.misc)))
})
test_that("The misc parameter won't impact gbm fitting when not passing censoring indicator - Cox Partial Hazards model", {
  require(survival)
  
  ## GIVEN
  # create some data
  N <- 3000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  f <- 0.5*sin(3*X1 + 5*X2^2 + mu/10)
  tt.surv <- rexp(N,exp(f))
  tt.cens <- rexp(N,0.5)
  delta <- as.numeric(tt.surv <= tt.cens)
  tt <- apply(cbind(tt.surv,tt.cens),1,min)
  
  # throw in some missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  # random weights if you want to experiment with them
  w <- rep(1,N)
  
  data <- data.frame(y=Surv(tt,delta),X1=X1,X2=X2,X3=X3)
  
  ## WHEN
  # Fit gbm with and without misc
  admisc <- runif(N)
  set.seed(15) 
  gbm.no.misc <- gbm.fit(x=data[,-1], y=data[,1], distribution = "coxph")
  
  set.seed(15)
  gbm.with.misc <- gbm.fit(x=data[,-1], y=data[,1], misc=admisc, distribution = "coxph")
  
  ## THEN
  expect_equal(gbm.no.misc, gbm.with.misc)
})
test_that("The misc parameter does not impact gbm fitting with - Gamma distribution", {
  # Create data
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  X4 <- ordered(sample(letters[1:6],N,replace=T))
  X5 <- factor(sample(letters[1:3],N,replace=T))
  X6 <- 3*runif(N)
  mu <- c(0,0,1,2)[as.numeric(X3)]
  
  SNR <- 10 # signal-to-noise ratio
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  sigma <- sqrt(var(Y)/SNR)
  Y <- Y + abs(rnorm(N,0,sigma))
  
  # create a bunch of missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  w <- rep(1,N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # Fit gbm with and without misc
  admisc <- runif(N)
  set.seed(15) 
  gbm.no.misc <- gbm.fit(x=data[,-1], y=data[,1], distribution = "gamma")
  
  set.seed(15)
  gbm.with.misc <- gbm.fit(x=data[,-1], y=data[,1], misc=admisc, distribution = "gamma")
  expect_equal(gbm.no.misc, gbm.with.misc)
})
test_that("The misc parameter does not impact gbm fitting with - Gaussian distribution", {
  # Create data
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  X4 <- ordered(sample(letters[1:6],N,replace=T))
  X5 <- factor(sample(letters[1:3],N,replace=T))
  X6 <- 3*runif(N)
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  SNR <- 10 # signal-to-noise ratio
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  sigma <- sqrt(var(Y)/SNR)
  Y <- Y + rnorm(N,0,sigma)
  
  # create a bunch of missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  w <- rep(1,N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # Fit gbm with and without misc
  admisc <- runif(N)
  set.seed(15) 
  gbm.no.misc <- gbm.fit(x=data[,-1], y=data[,1], distribution = "gaussian")
  
  set.seed(15)
  gbm.with.misc <- gbm.fit(x=data[,-1], y=data[,1], misc=admisc, distribution = "gaussian")
  expect_equal(gbm.no.misc, gbm.with.misc)
})
test_that("The misc parameter does not impact gbm fitting with - Huberized Hinge Loss", {
  # create some data
  set.seed(1)
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  
  # random weights if you want to experiment with them
  w <- rexp(N)
  w <- N*w/sum(w)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  
  # Fit gbm with and without misc
  admisc <- runif(N)
  set.seed(15) 
  gbm.no.misc <- gbm.fit(x=data[,-1], y=data[,1], distribution = "huberized")
  
  set.seed(15)
  gbm.with.misc <- gbm.fit(x=data[,-1], y=data[,1], misc=admisc, distribution = "huberized")
  
})
test_that("The misc parameter does not impact gbm fitting with - Laplace distribution", {
  # Create data
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  X4 <- ordered(sample(letters[1:6],N,replace=T))
  X5 <- factor(sample(letters[1:3],N,replace=T))
  X6 <- 3*runif(N)
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  SNR <- 10 # signal-to-noise ratio
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  sigma <- sqrt(var(Y)/SNR)
  Y <- Y + rnorm(N,0,sigma)
  
  # create a bunch of missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  w <- rep(1,N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # Fit gbm with and without misc
  admisc <- runif(N)
  set.seed(15) 
  gbm.no.misc <- gbm.fit(x=data[,-1], y=data[,1], distribution = "laplace")
  
  set.seed(15)
  gbm.with.misc <- gbm.fit(x=data[,-1], y=data[,1], misc=admisc, distribution = "laplace")
  expect_equal(gbm.no.misc, gbm.with.misc)
})
test_that("The misc parameter DOES impact gbm fitting with - Pairwise mrr & map", {
  # Create data
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
  Y.norm <- round(Y/(max(Y) + 0.001))
  
  data.train.mrr <- data.frame(Y.norm, query=query, X1, X2, X3)
  
  # Fit gbm with two misc's
  admisc.1 <- runif(N)
  admisc.2 <- runif(N)
  set.seed(15) 
  gbm.with.misc.mrr <- gbm.fit(x=data.train.mrr[,3:5], y=data.train.mrr[,1], misc=admisc.1, distribution=list(name='pairwise', metric="mrr", group='query'))
  gbm.with.misc.map <- gbm.fit(x=data.train.mrr[,3:5], y=data.train.mrr[,1], misc=admisc.1, distribution=list(name='pairwise', metric="map", group='query'))

  set.seed(15)
  gbm.with.misc.mrr.two <- gbm.fit(x=data.train.mrr[,3:5], y=data.train.mrr[,1], misc=admisc.2, distribution=list(name='pairwise', metric="mrr", group='query'))
  gbm.with.misc.map.two <- gbm.fit(x=data.train.mrr[,3:5], y=data.train.mrr[,1], misc=admisc.2, distribution=list(name='pairwise', metric="map", group='query'))
  expect_false(isTRUE(all.equal(gbm.with.misc.mrr, gbm.with.misc.mrr.two)))
  expect_false(isTRUE(all.equal(gbm.with.misc.map, gbm.with.misc.map.two)))
})
test_that("The misc parameter does not impact gbm fitting with - Poisson", {
  # Create data
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  X4 <- ordered(sample(letters[1:6],N,replace=T))
  X5 <- factor(sample(letters[1:3],N,replace=T))
  X6 <- 3*runif(N)
  mu <- c(0,0,1,2)[as.numeric(X3)]
  
  SNR <- 10 # signal-to-noise ratio
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  sigma <- sqrt(var(Y)/SNR)
  Y <- round(Y + abs(rnorm(N,0,sigma)))
  
  # create a bunch of missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  w <- rep(1,N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # Fit gbm with and without misc
  admisc <- runif(N)
  set.seed(15) 
  gbm.no.misc <- gbm.fit(x=data[,-1], y=data[,1], distribution = "poisson")
  
  set.seed(15)
  gbm.with.misc <- gbm.fit(x=data[,-1], y=data[,1], misc=admisc, distribution = "poisson")
  expect_equal(gbm.no.misc, gbm.with.misc)
})
test_that("The misc parameter does not impact gbm fitting with - Quantile Regressions", {
  # Create data
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  X4 <- ordered(sample(letters[1:6],N,replace=T))
  X5 <- factor(sample(letters[1:3],N,replace=T))
  X6 <- 3*runif(N)
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  SNR <- 10 # signal-to-noise ratio
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  sigma <- sqrt(var(Y)/SNR)
  Y <- Y + rnorm(N,0,sigma)
  
  # create a bunch of missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  w <- rep(1,N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # Fit gbm with and without misc
  admisc <- runif(N)
  set.seed(15) 
  gbm.no.misc <- gbm.fit(x=data[,-1], y=data[,1], distribution = list(name="quantile",alpha=0.95))
  
  set.seed(15)
  gbm.with.misc <- gbm.fit(x=data[,-1], y=data[,1], misc=admisc, distribution = list(name="quantile",alpha=0.95))
  expect_equal(gbm.no.misc, gbm.with.misc)
})
test_that("The misc parameter does not impact gbm fitting with - T distribution", {
  # Create data
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  X4 <- ordered(sample(letters[1:6],N,replace=T))
  X5 <- factor(sample(letters[1:3],N,replace=T))
  X6 <- 3*runif(N)
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  SNR <- 10 # signal-to-noise ratio
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  sigma <- sqrt(var(Y)/SNR)
  Y <- Y + rnorm(N,0,sigma)
  
  # create a bunch of missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  w <- rep(1,N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # Fit gbm with and without misc
  admisc <- runif(N)
  set.seed(15) 
  gbm.no.misc <- gbm.fit(x=data[,-1], y=data[,1], distribution = "tdist")
  
  set.seed(15)
  gbm.with.misc <- gbm.fit(x=data[,-1], y=data[,1], misc=admisc, distribution = "tdist")
  expect_equal(gbm.no.misc, gbm.with.misc)
})
test_that("The misc parameter does not impact gbm fitting with - Tweedie distribution", {
  # Create data
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  X4 <- ordered(sample(letters[1:6],N,replace=T))
  X5 <- factor(sample(letters[1:3],N,replace=T))
  X6 <- 3*runif(N)
  mu <- c(0,0,1,2)[as.numeric(X3)]
  
  SNR <- 10 # signal-to-noise ratio
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  sigma <- sqrt(var(Y)/SNR)
  Y <- Y + abs(rnorm(N,0,sigma))
  
  # create a bunch of missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  w <- rep(1,N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # Fit gbm with and without misc
  admisc <- runif(N)
  set.seed(15) 
  gbm.no.misc <- gbm.fit(x=data[,-1], y=data[,1], distribution = "tweedie")
  
  set.seed(15)
  gbm.with.misc <- gbm.fit(x=data[,-1], y=data[,1], misc=admisc, distribution = "tweedie")
  expect_equal(gbm.no.misc, gbm.with.misc)
})