####################
# Author: James Hickey
#
# Series of tests to check the calculation of the variable's 
# interaction strength
#
####################

context("Testing help functions - Interact")
test_that("check_and_set_num_trees throws a warning if given a vector", {
  # Given a vector for num_trees and a fake fit object
  gbm_fake_fit <- list()
  gbm_fake_fit$params$num_trees <- 2000
  num_trees <- c(200, 500)
  
  # When checking and setting the number of trees
  # Then a warning is thrown
  expect_warning(check_and_set_num_trees(gbm_fake_fit, num_trees), 
                 "length num_trees > 1: using first element")
})
test_that("check_and_set_num_trees throws an error if not given natural number for num_trees", {
  # Given a fake fit object and a num_trees which isn't a positive integer
  gbm_fake_fit <- list()
  gbm_fake_fit$params$num_trees <- 2000
  num_trees <- "wrong"
  
  # When checking and setting the number of trees
  # Then an error is thrown
  expect_error(check_and_set_num_trees(gbm_fake_fit, num_trees))
})
test_that("check_and_set_num_trees throws a warning if number of trees exceeds those in original fit", {
  # Given a fake fit object and num_trees > those fitted
  gbm_fake_fit <- list()
  gbm_fake_fit$params$num_trees <- 2000
  num_trees <- 5000
  
  # When checking and setting the number of trees
  # Then a warning is thrown
  expect_warning(check_and_set_num_trees(gbm_fake_fit, num_trees), 
                 paste("num_trees exceeds the number of trees in the model, ",
                       gbm_fake_fit$params$num_trees,". Using ", gbm_fake_fit$params$num_trees, " trees.", 
                       sep = ""))
  
})
test_that("check_and_set_num_trees returns 1st element of vector is < those in original fit", {
  # Given a vector for num_trees and a fake fit object
  gbm_fake_fit <- list()
  gbm_fake_fit$params$num_trees <- 2000
  num_trees <- c(200, 500)
  
  # When checking and setting the number of trees
  # Then returns first element of vector
  expect_equal(check_and_set_num_trees(gbm_fake_fit, num_trees), 
                num_trees[1])
})
test_that("check_and_set_num_trees returns those in original fit if number specified exceeds total in fit", {
  # Given a fake fit object and num_trees > those fitted
  gbm_fake_fit <- list()
  gbm_fake_fit$params$num_trees <- 2000
  num_trees <- 5000
  
  # When checking and setting the number of trees
  # Then returns number of trees in fit
  expect_equal(check_and_set_num_trees(gbm_fake_fit, num_trees), gbm_fake_fit$params$num_trees)
})
test_that("check_and_set_num_trees returns argument num_trees if it is a positive integer <= the total num_trees in original fit", {
  # Given a fake fit object and num_trees < those fitted
  gbm_fake_fit <- list()
  gbm_fake_fit$params$num_trees <- 2000
  num_trees <- 100
  
  # When checking and setting the number of trees
  # Then returns the num_trees passed in
  expect_equal(check_and_set_num_trees(gbm_fake_fit, num_trees), num_trees)
})
test_that("check_and_set_variables_indices throws an error if all strings but not all in fit's variable names", {
  # Given a fake fit object and a character vector with unrecognised variables
  gbm_fake_fit <- list()
  gbm_fake_fit$variables$var_names <- c("Var1", "Var2", "Var3")
  var_indices <- "Var4"
  
  # When checking and setting variable indices
  # Then an error is thrown
  expect_error(check_and_set_variables_indices(gbm_fake_fit, var_indices))
})
test_that("check_and_set_variables_indices throws a warning if indices not between 1 and the number of variables in fit", {
  # Given a fake fit object and index outside range
  gbm_fake_fit <- list()
  gbm_fake_fit$variables$var_names <- c("Var1", "Var2", "Var3")
  var_indices <- 4
  
  # When checking and setting variable indices
  # Then a warning will be thrown
  expect_warning(check_and_set_variables_indices(gbm_fake_fit, var_indices),
                 paste0("variables_indices must be between 1 and ", 
                        length(gbm_fake_fit$variables$var_names)))
})
test_that("check_and_set_variables_indices converts character vector input to integers", {
  # Given a fake fit object and variable indices that are strings
  gbm_fake_fit <- list()
  gbm_fake_fit$variables$var_names <- c("Var1", "Var2", "Var3")
  var_indices <- c("Var1", "Var2")
  
  # When checking and setting variables indices
  # Then the indices are converted to appropriate integers
  expect_equal(check_and_set_variables_indices(gbm_fake_fit, var_indices),
               c(1, 2))
})
test_that("check_and_set_variables_indices returns original vector if positive integers", {
  # Given a fake fit object and correct indices
  gbm_fake_fit <- list()
  gbm_fake_fit$variables$var_names <- c("Var1", "Var2", "Var3")
  var_indices <- c(2, 3)
  
  # When checking and setting variables indices
  # Then the original indices are returned
  expect_equal(check_and_set_variables_indices(gbm_fake_fit, var_indices),
               var_indices)
})
test_that("table_of_unique_values works correctly", {
  # Given data and correct variable indices
  N <- 100
  data <- data.frame("Var1"=rnorm(N), "Var2"=sample(1:5, N, replace=TRUE),
                     "Var3"=rnorm(N))
  var_indices <- c(1, 2)
  
  # When table_of_unique_values is called
  output <- table_of_unique_values(data, var_indices)
  
  # Then correct output produced
  corr_out <- unique(data[, var_indices, drop=FALSE])
  corr_out$num_levels_factors <- table(factor(apply(data[, var_indices,drop=FALSE],1,paste,collapse="\r"),
                                              levels=apply(corr_out, 1,paste,collapse="\r")))
  expect_equal(output, corr_out)
})
test_that("output of compute_preds_for_all_var_combinations is a list with the correct fields", {
  # Given data, a fit object, variable indices and all combinations of variables
  var_indices <- c(1, 2, 3)
  all_combinations_vars <- apply(expand.grid(rep(list(c(FALSE,TRUE)), length(var_indices)))[-1,], 1,
                                 function(x) as.numeric(which(x)))
  
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
  
  # random weights if you want to experiment with them
  w <- rep(1,N)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # fit initial model
  gbm_fit_obj <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
                     data=data,                   # dataset
                     var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                     distribution="Gaussian",     # bernoulli, adaboost, gaussian, poisson, coxph, or
                     # list(name="quantile",alpha=0.05) for quantile regression
                     n.trees=2000,                 # number of trees
                     shrinkage=0.005,             # shrinkage or learning rate, 0.001 to 0.1 usually work
                     interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc
                     bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
                     train.fraction = 0.5,        # fraction of data for training, first train.fraction*N used for training
                     mFeatures = 3,               # Number of features to consider at each node.
                     n.minobsinnode = 10,         # minimum number of obs needed in each node
                     keep.data=TRUE,
                     cv.folds=10,                 # do 10-fold cross-validation
                     verbose = FALSE)             # don't print progress
  
  # When calculating predictions for all combinations of variables
  # Then has correct fields
  expect_equal(names(compute_preds_for_all_var_combinations(data, gbm_fit_obj,
                                                      all_combinations_vars, var_indices, 1000)[[1]]),
  c("data", "num_levels_factors", "preds", "sign"))
})
test_that("output of compute_pred_for_all_var_combinations is of the correct size", {
  # Given data, a fit object, variable indices and all combinations of variables
  var_indices <- c(1, 2, 3)
  all_combinations_vars <- apply(expand.grid(rep(list(c(FALSE,TRUE)), length(var_indices)))[-1,], 1,
                                 function(x) as.numeric(which(x)))
  
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
  
  # random weights if you want to experiment with them
  w <- rep(1,N)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # fit initial model
  gbm_fit_obj <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
                     data=data,                   # dataset
                     var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                     distribution="Gaussian",     # bernoulli, adaboost, gaussian, poisson, coxph, or
                     # list(name="quantile",alpha=0.05) for quantile regression
                     n.trees=2000,                 # number of trees
                     shrinkage=0.005,             # shrinkage or learning rate, 0.001 to 0.1 usually work
                     interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc
                     bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
                     train.fraction = 0.5,        # fraction of data for training, first train.fraction*N used for training
                     mFeatures = 3,               # Number of features to consider at each node.
                     n.minobsinnode = 10,         # minimum number of obs needed in each node
                     keep.data=TRUE,
                     cv.folds=10,                 # do 10-fold cross-validation
                     verbose = FALSE)             # don't print progress
  
  # When calculating predictions for all combinations of variables
  # Then has correct size
  expect_equal(length(compute_preds_for_all_var_combinations(data, 
                                                             gbm_fit_obj,
                                                             all_combinations_vars, 
                                                             var_indices, 
                                                             num_trees = 1000)),
               length(all_combinations_vars))
})
test_that("sign field of output of compute_pred_for_all_var_combinations consists of +/-1s", {
  # Given data, a fit object, variable indices and all combinations of variables
  var_indices <- c(1, 2, 3)
  all_combinations_vars <- apply(expand.grid(rep(list(c(FALSE,TRUE)), length(var_indices)))[-1,], 1,
                                 function(x) as.numeric(which(x)))
  
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
  
  # random weights if you want to experiment with them
  w <- rep(1,N)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # fit initial model
  gbm_fit_obj <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
                     data=data,                   # dataset
                     var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                     distribution="Gaussian",     # bernoulli, adaboost, gaussian, poisson, coxph, or
                     # list(name="quantile",alpha=0.05) for quantile regression
                     n.trees=2000,                 # number of trees
                     shrinkage=0.005,             # shrinkage or learning rate, 0.001 to 0.1 usually work
                     interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc
                     bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
                     train.fraction = 0.5,        # fraction of data for training, first train.fraction*N used for training
                     mFeatures = 3,               # Number of features to consider at each node.
                     n.minobsinnode = 10,         # minimum number of obs needed in each node
                     keep.data=TRUE,
                     cv.folds=10,                 # do 10-fold cross-validation
                     verbose = FALSE)             # don't print progress
  
  # When calculating predictions for all combinations of variables
  preds_for_combs <- compute_preds_for_all_var_combinations(data, 
                                                            gbm_fit_obj,
                                                            all_combinations_vars, 
                                                            var_indices, 
                                                            num_trees = 1000)
  # Then sign field is +-1s
  expect_true(all(preds_for_combs[[1]]$sign %in% c(1, -1)))
})

context("testing interact.GBMFit")
test_that("if data not a data.frame or matrix then an error is thrown", {
  # Given data and a fit object 
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
  
  # random weights if you want to experiment with them
  w <- rep(1,N)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # fit initial model
  gbm_fit_obj <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
                     data=data,                   # dataset
                     var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                     distribution="Gaussian",     # bernoulli, adaboost, gaussian, poisson, coxph, or
                     # list(name="quantile",alpha=0.05) for quantile regression
                     n.trees=2000,                 # number of trees
                     shrinkage=0.005,             # shrinkage or learning rate, 0.001 to 0.1 usually work
                     interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc
                     bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
                     train.fraction = 0.5,        # fraction of data for training, first train.fraction*N used for training
                     mFeatures = 3,               # Number of features to consider at each node.
                     n.minobsinnode = 10,         # minimum number of obs needed in each node
                     keep.data=TRUE,
                     cv.folds=10,                 # do 10-fold cross-validation
                     verbose = FALSE)             # don't print progress
  
  # When calculating the model interactions with data
  # not in data frame or matrix
  # Then an error is thrown
  expect_error(interact(gbm_fit_obj, as.list(data)), 
               "data argument should be a data.frame or matrix")
})
test_that("if var_indices not a vector of characters or integers then an error is thrown", {
  # Given data and a fit object 
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
  
  # random weights if you want to experiment with them
  w <- rep(1,N)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # fit initial model
  gbm_fit_obj <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
                     data=data,                   # dataset
                     var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                     distribution="Gaussian",     # bernoulli, adaboost, gaussian, poisson, coxph, or
                     # list(name="quantile",alpha=0.05) for quantile regression
                     n.trees=2000,                 # number of trees
                     shrinkage=0.005,             # shrinkage or learning rate, 0.001 to 0.1 usually work
                     interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc
                     bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
                     train.fraction = 0.5,        # fraction of data for training, first train.fraction*N used for training
                     mFeatures = 3,               # Number of features to consider at each node.
                     n.minobsinnode = 10,         # minimum number of obs needed in each node
                     keep.data=TRUE,
                     cv.folds=10,                 # do 10-fold cross-validation
                     verbose = FALSE)             # don't print progress
  
  # When calculating interactions using var_indices that aren't a
  # vector of integers or character vectors
  var_ind1 <- Inf
  var_ind2 <- NA
  var_ind3 <- c(1, "Hello")
  
  # Then error is thrown
  expect_error(interact(gbm_fit_obj, data, var_ind1), 
               "Variables indices must be a vector of integers or characters")
  expect_error(interact(gbm_fit_obj, data, var_ind2), 
               "Variables indices must be a vector of integers or characters")
  expect_error(interact(gbm_fit_obj, data, var_ind3), 
               "Variables indices must be a vector of integers or characters")
})
test_that("Error thrown if length var_indices exceeds interaction depth of fit", {
  # Given data and a fit object 
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
  
  # random weights if you want to experiment with them
  w <- rep(1,N)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # fit initial model
  gbm_fit_obj <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
                     data=data,                   # dataset
                     var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                     distribution="Gaussian",     # bernoulli, adaboost, gaussian, poisson, coxph, or
                     # list(name="quantile",alpha=0.05) for quantile regression
                     n.trees=2000,                 # number of trees
                     shrinkage=0.005,             # shrinkage or learning rate, 0.001 to 0.1 usually work
                     interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc
                     bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
                     train.fraction = 0.5,        # fraction of data for training, first train.fraction*N used for training
                     mFeatures = 3,               # Number of features to consider at each node.
                     n.minobsinnode = 10,         # minimum number of obs needed in each node
                     keep.data=TRUE,
                     cv.folds=10,                 # do 10-fold cross-validation
                     verbose = FALSE)             # don't print progress
  
  # When evaluating the interactions for more than variables than interaction depth
  # Then an error is thrown
  expect_error(interact(gbm_fit_obj, data, 1:4))
})
test_that("Given correct inputs will return H-statistic < 1", {
  # Given data and a fit object 
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
  
  # random weights if you want to experiment with them
  w <- rep(1,N)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # fit initial model
  gbm_fit_obj <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
                    data=data,                   # dataset
              var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
              distribution="Gaussian",     # bernoulli, adaboost, gaussian, poisson, coxph, or
              # list(name="quantile",alpha=0.05) for quantile regression
              n.trees=2000,                 # number of trees
              shrinkage=0.005,             # shrinkage or learning rate, 0.001 to 0.1 usually work
              interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc
              bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
              train.fraction = 0.5,        # fraction of data for training, first train.fraction*N used for training
              mFeatures = 3,               # Number of features to consider at each node.
              n.minobsinnode = 10,         # minimum number of obs needed in each node
              keep.data=TRUE,
              cv.folds=10,                 # do 10-fold cross-validation
              verbose = FALSE)             # don't print progress
  
  # When calculating the interactions 
  int_stats <- interact(gbm_fit_obj, data=data, var_indices=1:2)
  
  # Then H-statistic is correct and < 1
  expect_true(is.double(int_stats))
  expect_true(int_stats < 1)
})