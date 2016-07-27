####################
# Author: James Hickey
#
# Series of tests for perf_pairwise functionality
#
####################

context("Testing perf_pairwise")
test_that("Metric ndcg", {
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
  
  fit <- gbm2(Y~X1+X2+X3,          # formula
             data=data,     # dataset
             distribution=dist,
             train_params=params,
             keep_gbm_data=TRUE,      # store copy of input data in model
             cv_folds=5,          # number of cross validation folds
             is_verbose = FALSE ,    # don't print progress
             par_details=gbmParallel())  
  
  
})