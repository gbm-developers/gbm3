####################
# Author: James Hickey
#
# Series of tests to check other, less used, distributions can be run to 
# completion
#
####################

context("Testing runs of other distributions")

test_that("Can fit GBM with - AdaBoost", {
  # Given appropriate data
  set.seed(1)
  
  # create some data
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
  offset <- rep(0, N)
  
  # Set up for new API
  params <- training_params(num_trees=20, interaction_depth=3,
                            min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5,
                            id=seq(nrow(data)), num_train=N/2,
                            num_features=3)
  
  # When GBM fit with AdaBoost
  dist <- gbm_dist("AdaBoost")
  
  # Then no error is thrown
  expect_error(gbmt(Y~X1+X2+X3, data=data,
                    distribution=dist, weights=w, offset=offset,
                    train_params=params, var_monotone=c(0, 0, 0), 
                    keep_gbm_data=TRUE, cv_folds=5, is_verbose = FALSE),
               NA)
})

test_that("Can fit GBM with - Gamma", {
  # Given appropriate data
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
  Y <- abs(Y + rnorm(N,0,sigma))
  
  # create a bunch of missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  w <- rep(1,N)
  offset <- rep(0, N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  
  # Set up for new API
  params <- training_params(num_trees=20, interaction_depth=3,
                            min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5,
                            id=seq(nrow(data)), num_train=N/2,
                            num_features=6)
  # When GBM fit with Gamma
  dist <- gbm_dist("Gamma")
  
  # Then no error is thrown
  expect_error(gbmt(Y~X1+X2+X3+X4+X5+X6, data=data,
                    distribution=dist, weights=w, offset=offset,
                    train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), 
                    keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE),
               NA)
})

test_that("Can fit GBM with - Huberized", {
  # Given appropriate data
  set.seed(1)
  
  # create some data
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
  offset <- rep(0, N)
  
  # Set up for new API
  params <- training_params(num_trees=20, interaction_depth=3,
                            min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5,
                            id=seq(nrow(data)), num_train=N/2,
                            num_features=3)
  
  # When GBM fit with Huberized
  dist <- gbm_dist("Huberized")
  
  # Then no error is thrown
  expect_error(gbmt(Y~X1+X2+X3, data=data, distribution=dist,
                    weights=w, offset=offset,
                    train_params=params, var_monotone=c(0, 0, 0), 
                    keep_gbm_data=TRUE, cv_folds=5, is_verbose = FALSE),
               NA)
})

test_that("Can fit GBM with - Laplace", {
  # Given appropriate data
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
  params <- training_params(num_trees=20, interaction_depth=3,
                            min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5,
                            id=seq(nrow(data)), num_train=N/2,
                            num_features=6)
  # When GBM fit with Laplace
  dist <- gbm_dist("Laplace")
  
  # Then no error is thrown
  expect_error(gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist,
                    weights=w, offset=offset,
                    train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), 
                    keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE),
               NA)
})

test_that("Can fit GBM with - Pairwise", {
  skip("Skipping pairwise")
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
  Y_M <- Y
  Y_M[Y_M >= 1] <- 1
  Y_M[Y_M < 1] <- 0
  
  data <- data.frame(Y, query=query, X1, X2, X3)
  data_M <- data.frame(Y_M, query=query, X1, X2, X3)
  
  params <- training_params(num_trees = 20, num_train = nrow(data),
                            id=seq_len(nrow(data)),
                            interaction_depth = 3)
  
  # When fitting all Pairwise distributions
  dist <- gbm_dist("Pairwise", metric="ndcg", group="query")
  dist_2 <- gbm_dist("Pairwise", metric="conc", group="query")
  dist_3 <- gbm_dist("Pairwise", metric="mrr", group="query")
  dist_4 <- gbm_dist("Pairwise", metric="map", group="query")

  # Then runs without error
  expect_error(gbmt(Y~X1+X2+X3,          
              data=data,    
              distribution=dist,
              train_params=params,
              keep_gbm_data=TRUE,     
              cv_folds=5,   
              is_verbose = FALSE ,    
              par_details=gbmParallel()), NA)
  expect_error(gbmt(Y~X1+X2+X3,          
                    data=data,    
                    distribution=dist_2,
                    train_params=params,
                    keep_gbm_data=TRUE,     
                    cv_folds=5,   
                    is_verbose = FALSE ,    
                    par_details=gbmParallel()), NA)
  expect_error(gbmt(Y_M~X1+X2+X3,          
                    data=data_M,    
                    distribution=dist_3,
                    train_params=params,
                    keep_gbm_data=TRUE,     
                    cv_folds=5,   
                    is_verbose = FALSE ,    
                    par_details=gbmParallel()), NA)
  expect_error(gbmt(Y_M~X1+X2+X3,          
                    data=data_M,    
                    distribution=dist_4,
                    train_params=params,
                    keep_gbm_data=TRUE,     
                    cv_folds=5,   
                    is_verbose = FALSE ,    
                    par_details=gbmParallel()), NA)
  
})

test_that("Can fit GBM with - Poisson", {
  # Given appropriate data
  set.seed(1)
  
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p) + 1
  
  # random weights if you want to experiment with them
  w <- rexp(N)
  w <- N*w/sum(w)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  offset <- rep(0, N)
  
  # Set up for new API
  params <- training_params(num_trees=20, interaction_depth=3,
                            min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5,
                            id=seq(nrow(data)), num_train=N/2, num_features=3)
  
  # When GBM fit with Poisson
  dist <- gbm_dist("Poisson")
  
  # Then no error is thrown
  expect_error(gbmt(Y~X1+X2+X3, data=data, distribution=dist,
                    weights=w, offset=offset,
                    train_params=params, var_monotone=c(0, 0, 0), 
                    keep_gbm_data=TRUE, cv_folds=5, is_verbose = FALSE),
               NA)
})

test_that("Can fit GBM with - Quantile", {
  ## Given appropriate data
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
  params <- training_params(num_trees=20, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  # When GBM fit with Quantile
  dist <- gbm_dist("Quantile")
  
  # Then no error is thrown
  expect_error(gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist,
                    weights=w, offset=offset,
                    train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), 
                    keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE),
               NA)
})

test_that("Can fit GBM with - TDist", {
  # Given appropriate data
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
  params <- training_params(num_trees=20, interaction_depth=3,
                            min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5,
                            id=seq(nrow(data)), num_train=N/2,
                            num_features=6)
  # When GBM fit with TDist
  dist <- gbm_dist("TDist")
  
  # Then no error is thrown
  expect_error(gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist,
                    weights=w, offset=offset,
                    train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), 
                    keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE),
               NA)
})

test_that("Can fit GBM with - Tweedie", {
  # Given appropriate data
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
  Y <- abs(Y + rnorm(N,0,sigma))
  
  # create a bunch of missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  w <- rep(1,N)
  offset <- rep(0, N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  
  # Set up for new API
  params <- training_params(num_trees=20, interaction_depth=3,
                            min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5,
                            id=seq(nrow(data)), num_train=N/2,
                            num_features=6)
  # When GBM fit with Gamma
  dist <- gbm_dist("Tweedie")
  
  # Then no error is thrown
  expect_error(gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist,
                    weights=w, offset=offset,
                    train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), 
                    keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE),
               NA)
})
