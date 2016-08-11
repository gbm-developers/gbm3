####################
# Author: James Hickey
#
# Series of test to validate conversion to old gbm class
#
####################

context("Testing conversion of GBMFit to gbm object")
test_that("to_old_gbm produces an object with the right class", {
  set.seed(1)
  
  # create some data
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
  offset <- rep(0, N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  
  # Set up for new API
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE)
  # When converted to old form
  # Then has right class
  expect_equal(class(to_old_gbm(fit)), "gbm")
})
test_that("to_old_gbm output has correct fields", {
  set.seed(1)
  
  # create some data
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
  offset <- rep(0, N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  
  # Set up for new API
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE)
  # When converted to old form
  # Then has right fields
  expect_equal(names(to_old_gbm(fit)), c("initF", "train.error", "valid.error", "trees",
                                         "c.splits", "oobag.improve", "fit", "bag.fraction",
                                         "distribution", "interaction.depth", "n.minobsinnode",
                                         "n.trees", "nTrain", "nTrainPats", "patient.id",
                                         "mFeatures", "train.fraction", "response.name", "shrinkage",
                                         "var.levels", "var.monotone", "var.names", "var.type", "verbose",
                                         "strata", "sorted", "prior.node.coeff.var", "data", "cv.folds", "cv.error", "Terms",
                                         "call", "m"))
})
test_that("to_old_gbm gives correct message on calling", {
  set.seed(1)
  
  # create some data
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
  offset <- rep(0, N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  
  # Set up for new API
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE)
  # When converted to old form
  # Then correct message given
  expect_message(to_old_gbm(fit), "Converted to old gbm object - this will not work with new packages functions")
})
test_that("to_old_gbm has data if keep_gbm_data was specified and it is correct", {
  set.seed(1)
  
  # create some data
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
  offset <- rep(0, N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  
  # Set up for new API
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE)
  # When converted to old form
  # Then data stored is correct
  expect_equal(to_old_gbm(fit)$data$y, fit$gbm_data_obj$y)
  expect_equal(to_old_gbm(fit)$data$x, fit$gbm_data_obj$x)
  expect_equal(to_old_gbm(fit)$data$w, fit$gbm_data_obj$w)
  expect_equal(to_old_gbm(fit)$data$offset, fit$gbm_data_obj$offset)
  expect_true(is.na(to_old_gbm(fit)$data$Misc))
})
test_that("to_old_gbm has correct default prior.node.coeff.var if not CoxPH", {
  set.seed(1)
  
  # create some data
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
  offset <- rep(0, N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  
  # Set up for new API
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE)
  # When to_old_gbm called
  # prior.node.coeff.var has correct default value
  expect_equal(to_old_gbm(fit)$prior.node.coeff.var, 1000)
})
test_that("to_old_gbm has correct data and time ordering for CoxPH fit", {
  # Require Surv to be available
  require(survival)
  
  # create some data
  set.seed(1)
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
  data <- data.frame(tt=tt,delta=delta,X1=X1,X2=X2,X3=X3)
  
  # Put into new API
  dist <- gbm_dist("CoxPH", prior_node_coeff_var=10)
  params <- training_params(num_trees=3000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=3)
  
  # fit initial model
  gbm1 <- gbmt(Surv(tt,delta)~X1+X2+X3, data=data, distribution=dist, weights=w, offset=rep(0, N),
               train_params=params, var_monotone=c(0, 0, 0), keep_gbm_data=TRUE, cv_folds=5, is_verbose=FALSE)
  # When to old called
  # Then tied.times.method, prior.node.coeff.var and i.order is specified
  expect_equal(to_old_gbm(gbm1)$tied.times.method, "efron")
  expect_equal(to_old_gbm(gbm1)$data$i.order, gbm1$distribution$time_order)
  expect_equal(to_old_gbm(gbm1)$prior.node.coeff.var, 10)
})
test_that("to_old_gbm has group ordering for Pairwise fit", {
  # Given data and a fitted pairwise
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
  params <- training_params(num_trees = 2000, num_train = nrow(data), id=seq_len(nrow(data)),
                            interaction_depth = 3)
  
  fit <- gbmt(Y~X1+X2+X3,          # formula
              data=data,     # dataset
              distribution=dist,
              train_params=params,
              keep_gbm_data=TRUE,      # store copy of input data in model
              cv_folds=5,          # number of cross validation folds
              is_verbose = FALSE ,    # don't print progress
              par_details=gbmParallel())  
  
  # When fit converted to old form
  # Then it has ord.group field
  expect_equal(to_old_gbm(fit)$ord.group, fit$distribution$group_order)
})