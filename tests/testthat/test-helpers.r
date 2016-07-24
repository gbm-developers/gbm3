####################
# Author: James Hickey
#
# Series of tests to check the helper functions
# within the package
#
####################

context("Testing checks on GBM S3 objects")
test_that("check_if_gbm_data throws an error if given an object not of class GBMData", {
  # Given an object not of class GBMData
  obj <- list()
  
  # When passed to check_if_gbm_data
  # Then an error is thrown
  expect_error(check_if_gbm_data(obj))
})

test_that("check_if_gbm_data does not throw an error if given an object of class GBMData", {
  # Given an object of class GBMData
  obj <- list()
  class(obj) <- "GBMData"
  
  # When passed to check_if_gbm_data
  # Then no error is thrown
  expect_error(check_if_gbm_data(obj), NA)
})

test_that("check_if_gbm_dist throws an error if given an object not of class GBMDist", {
  # Given an object not of class GBMDist
  obj <- list()
  
  # When passed to check_if_gbm_dist
  # Then an error is thrown
  expect_error(check_if_gbm_dist(obj))
})

test_that("check_if_gbm_dist does not throw an error if given an object of class GBMDist", { 
  # Given an object of class GBMDist
  obj <- list()
  class(obj) <- "GBMDist"
  
  # When passed to check_if_gbm_dist
  # Then no error is thrown
  expect_error(check_if_gbm_dist(obj), NA)  
})

test_that("check_if_gbm_train_params throws an error if given object not of class GBMTrainParams", {
  # Given an object not of class GBMTrainParams
  obj <- list()
  
  # When passed to check_if_gbm_train_params
  # Then an error is thrown
  expect_error(check_if_gbm_train_params(obj))
})

test_that("check_if_gbm_train_params does not throw an error if given object of class GBMTrainParams", {
  # Given an object of class GBMTrainParams
  obj <- list()
  class(obj) <- "GBMTrainParams"
  
  # When passed to check_if_gbm_train_params
  # Then no error is thrown
  expect_error(check_if_gbm_train_params(obj), NA)
})

test_that("check_if_gbm_fit throws an error if given an object not of class GBMFit", {
  # Given an object not of class GBMFit
  obj <- list()
  
  # When passed to check_if_gbm_fit
  # Then an error is thrown
  expect_error(check_if_gbm_fit(obj))
})

test_that("check_if_gbm_fit does not throw an error if given an object of class GBMFit", {
  # Given an object of class GBMFit
  obj <- list()
  class(obj) <- "GBMFit"
  
  # When passed to check_if_gbm_fit
  # Then no error is thrown
  expect_error(check_if_gbm_fit(obj), NA)
})

test_that("check_if_gbm_var_container throws an error if given an object not of class GBMVarCont", {
  # Given an object not of class GBMVarCont
  obj <- list()
  
  # When passed to check_if_gbm_var_container
  # Then an error is thrown
  expect_error(check_if_gbm_var_container(obj))
})

test_that("check_if_gbm_var_container does not throw an error if given an object of class GBMVarCont", {
  # Given an object of class GBMVarCont
  obj <- list()
  class(obj) <- "GBMVarCont"
  
  # When passed to check_if_gbm_var_container
  # Then no error is thrown
  expect_error(check_if_gbm_var_container(obj), NA)
})


context("Testing checks on inputs for creation of S3 objects")

test_that("check_weights returns a vector of 1s if passed an obj of length 0", {
  # Given an empty vector
  w <- c()
  n <- 100
  
  # When check_weights is called 
  w <- check_weights(w, n)
  
  # Then returns a vector of 1s
  expect_equal(w, rep(1, n))
})

test_that("check_weights throws an error if any not a double", {
  # Given a vector of weights
  n <- 100
  w <- rep(1, 100)
  
  # When one weight is not a double
  w[1] <- NA
  
  # Then check_weights will throw an error
  expect_error(check_weights(w, n))
})

test_that("check_weights throws an error if any are < 0", {
  # Given a vector of weights
  w <- rep(1, 100)
  
  # When one weight is set < 0
  w[1] <- -1
  
  # Then check_weights will throw an error
  expect_error(check_weights(w, 100))
})

test_that("check_interaction_depth throws error if < 1 or > 49", {
  # Given interaction depths of 0 and 50
  id1 <- 0
  id2 <- 50
  
  # When they're checked
  # Then error
  expect_error(check_interaction_depth(id1))
  expect_error(check_interaction_depth(id2))
})

test_that("checkMissing throws error if NaN in predictors", {
  # Given predictors and response
  N <- 100
  x <- runif(N)
  y <- runif(N)
  
  # When x has a NaN in it
  x[1] <- NaN
  
  # Then error thrown by checkMissing
  expect_error(checkMissing(x, y))
})
test_that("checkMissing throws error if missing value in response", {
  # Given predictors and response
  N <- 100
  x <- runif(N)
  y <- runif(N)
  
  # When y has a missing value
  y[1] <- NA
  
  # Then error thrown by checkMissing
  expect_error(checkMissing(x, y))
})

test_that("warnNoVariation warns if a variable has no variation", {
  x <- c(1.2, 1.2)
  
  expect_warning(warnNoVariation(x, 1, 'test'),
                 "variable 1: test has no variation",
                 fixed=TRUE)
})

test_that("warnNoVariation passes OK if variable does vary", {
  x <- c(0, 1)
  
  expect_warning(warnNoVariation(x, 1, 'test'), regexp=NA)
})

test_that("warnNoVariation passes OK if variable does vary (with NA)", {
  x <- c(0, 1, NA)
  
  expect_warning(warnNoVariation(x, 1, 'test'), regexp=NA)
})

test_that("get_var_names gets colnames if passed a matrix", {
  # Given a matrix with 3 columns and names
  x <- matrix(ncol=3)
  colnames(x) <- c("V1", "V2", "V3")
  
  # When passed to get_var_names
  # Then column names returned
  expect_equal(get_var_names(x), colnames(x))
})

test_that("get_var_names gets the names if passed a data.frame", {
  # Given a data.frame with 3 columns and names
  x <- data.frame(X1=NA, X2=NA, X3=NA)
  
  # When passed to get_var_names
  # Then names returned
  expect_equal(get_var_names(x), names(x))
})

test_that("check_sanity only throws an error if nrows(x) is not equal to the length of y", {
  # Given two sets of predictors one of length y, the other not
  N <- 100
  x1 <- data.frame(runif(N))
  x2 <- data.frame(runif(N-1))
  y <- runif(N)
  
  # When check_sanity is called
  # Then error thrown for call where lengths are not equal
  expect_error(check_sanity(x1, y), NA)
  expect_error(check_sanity(x2, y))
})

test_that("check_if_natural_number throws error if not passed a whole number >=1" , {
  # When check_if_natural_number is not passed a whole number >= 1
  # Then an error is thrown
  expect_error(check_if_natural_number(-1, "Arg Name"))
  expect_error(check_if_natural_number(c(1, 2), "Arg Name"))
  expect_error(check_if_natural_number(1.4, "Arg Name"))
  expect_error(check_if_natural_number(Inf, "Arg Name"))
  expect_error(check_if_natural_number(list(), "Arg Name"))
  expect_error(check_if_natural_number(NA, "Arg Name"))
  expect_error(check_if_natural_number(NaN, "Arg Name"))
})

test_that("convertY changes factors with 2 levels to numeric", {
  # Given a vector of 2-level factors
  N <- 100
  y <- as.factor(sample(c(0, 1), N, replace=TRUE))
  
  # When convertY is called
  # Then converted to numeric
  expect_equal(convertY(y), as.numeric(y==levels(y)[2]))
})

test_that("convertY does nothing to factors with levels != 2", {
  # Given of a vector of 3-level factors
  N <- 100
  y <- as.factor(sample(c(0, 1, 2)), N, replace=TRUE)
  
  # When convertY is called
  # Then remains the same
  expect_equal(convertY(y), y)
})
test_that("check_var_type throws error when passed excessive levels in x", {
  
  testExcess <- data.frame(
    y = sample(c(0,1), 1025, replace = TRUE)
    ,x1 = runif(1025)
    ,x2 = factor(1:1025)
  )
  
  expect_error(
    check_var_type(x = testExcess[,c('x1', 'x2')]
             , y = testExcess$y)
    ,'gbm does not currently handle categorical variables with more than 1024 levels\\. Variable 2\\: x2 has 1025 levels\\.')
  
})

test_that("check_var_type throws an error when passed Inacceptable classes", {
  
  testClasses <- data.frame(
    y = sample(c(0,1), 15, replace = TRUE)
    ,x1 = runif(15)
    ,x2 = seq(as.Date('2015-01-01'), as.Date('2015-01-15'), 'days')
  )
  
  expect_error(
    check_var_type(x = testClasses[,c('x1', 'x2')]
             , y = testClasses$y)
    ,'variable 2\\: x2 is not of type - numeric, ordered or factor\\.')
  
})

test_that("check_offset default returns a vector of 0s when offset set to NULL - irrespective of distribution", {
  # Given an offset=NULL and responses (and all distributions)
  N < - 100
  y <- runif(N)
  offset <- NULL
  dist_1 <- gbm_dist("AdaBoost")
  dist_2 <- gbm_dist("Bernoulli")
  dist_3 <- gbm_dist("CoxPH")
  dist_4 <- gbm_dist("Gamma")
  dist_5 <- gbm_dist("Gaussian")
  dist_6 <- gbm_dist("Huberized")
  dist_7 <- gbm_dist("Laplace")
  dist_8 <- gbm_dist("Pairwise")
  dist_9 <- gbm_dist("Poisson")
  dist_10 <- gbm_dist("Quantile")
  dist_11 <- gbm_dist("TDist")
  dist_12 <- gbm_dist("Tweedie")
  
  # Then check_offset returns a vector of zeros
  # equal to length of response
  expect_error(check_offset(offset, y, dist_1), rep(0, length(y)))
  expect_error(check_offset(offset, y, dist_2), rep(0, length(y)))
  expect_error(check_offset(offset, y, dist_3), rep(0, length(y)))
  expect_error(check_offset(offset, y, dist_4), rep(0, length(y)))
  expect_error(check_offset(offset, y, dist_5), rep(0, length(y)))
  expect_error(check_offset(offset, y, dist_6), rep(0, length(y)))
  expect_error(check_offset(offset, y, dist_7), rep(0, length(y)))
  expect_error(check_offset(offset, y, dist_8), rep(0, length(y)))
  expect_error(check_offset(offset, y, dist_9), rep(0, length(y)))
  expect_error(check_offset(offset, y, dist_10), rep(0, length(y)))
  expect_error(check_offset(offset, y, dist_11), rep(0, length(y)))
  expect_error(check_offset(offset, y, dist_12), rep(0, length(y)))
})

test_that("check_offset throws an error length of offset does not equal the length of the response - and not CoxPH", {
  # Given an offset and vector of responses
  # offset is different length to responses
  N < - 100
  y <- runif(N)
  offset <- runif(N-2)
  dist_1 <- gbm_dist("AdaBoost")
  dist_2 <- gbm_dist("Bernoulli")
  dist_3 <- gbm_dist("CoxPH")
  dist_4 <- gbm_dist("Gamma")
  dist_5 <- gbm_dist("Gaussian")
  dist_6 <- gbm_dist("Huberized")
  dist_7 <- gbm_dist("Laplace")
  dist_8 <- gbm_dist("Pairwise")
  dist_9 <- gbm_dist("Poisson")
  dist_10 <- gbm_dist("Quantile")
  dist_11 <- gbm_dist("TDist")
  dist_12 <- gbm_dist("Tweedie")
  
  # Then check_offset throws an error - if not CoxPH
  expect_error(check_offset(offset, y, dist_1))
  expect_error(check_offset(offset, y, dist_2))
  expect_error(check_offset(offset, y, dist_3), NA)
  expect_error(check_offset(offset, y, dist_4))
  expect_error(check_offset(offset, y, dist_5))
  expect_error(check_offset(offset, y, dist_6))
  expect_error(check_offset(offset, y, dist_7))
  expect_error(check_offset(offset, y, dist_8))
  expect_error(check_offset(offset, y, dist_9))
  expect_error(check_offset(offset, y, dist_10))
  expect_error(check_offset(offset, y, dist_11))
  expect_error(check_offset(offset, y, dist_12))
})

test_that("check_offset throws an error if the offset contains a NA", {
  # Given an offset and vector of responses - irrespective of distribution
  N < - 100
  y <- runif(N)
  offset <- runif(N)
  dist_1 <- gbm_dist("AdaBoost")
  dist_2 <- gbm_dist("Bernoulli")
  dist_3 <- gbm_dist("CoxPH")
  dist_4 <- gbm_dist("Gamma")
  dist_5 <- gbm_dist("Gaussian")
  dist_6 <- gbm_dist("Huberized")
  dist_7 <- gbm_dist("Laplace")
  dist_8 <- gbm_dist("Pairwise")
  dist_9 <- gbm_dist("Poisson")
  dist_10 <- gbm_dist("Quantile")
  dist_11 <- gbm_dist("TDist")
  dist_12 <- gbm_dist("Tweedie")
  
  
  # When an elemenet of offset is NA
  offset[1] <- NA
  
  # Then check_offset throws an error
  expect_error(check_offset(offset, y, dist_1))
  expect_error(check_offset(offset, y, dist_2))
  expect_error(check_offset(offset, y, dist_3), NA)
  expect_error(check_offset(offset, y, dist_4))
  expect_error(check_offset(offset, y, dist_5))
  expect_error(check_offset(offset, y, dist_6))
  expect_error(check_offset(offset, y, dist_7))
  expect_error(check_offset(offset, y, dist_8))
  expect_error(check_offset(offset, y, dist_9))
  expect_error(check_offset(offset, y, dist_10))
  expect_error(check_offset(offset, y, dist_11))
  expect_error(check_offset(offset, y, dist_12))
})

test_that("check_offset throws an error if the offset contains a non-numeric", {
  # Given an offset and vector of responses - irrespective of distribution
  N < - 100
  y <- runif(N)
  offset <- runif(N)
  dist_1 <- gbm_dist("AdaBoost")
  dist_2 <- gbm_dist("Bernoulli")
  dist_3 <- gbm_dist("CoxPH")
  dist_4 <- gbm_dist("Gamma")
  dist_5 <- gbm_dist("Gaussian")
  dist_6 <- gbm_dist("Huberized")
  dist_7 <- gbm_dist("Laplace")
  dist_8 <- gbm_dist("Pairwise")
  dist_9 <- gbm_dist("Poisson")
  dist_10 <- gbm_dist("Quantile")
  dist_11 <- gbm_dist("TDist")
  dist_12 <- gbm_dist("Tweedie")
  
  
  # When an elemenet of offset is NA
  offset[1] <- TRUE
  
  # Then check_offset throws an error
  expect_error(check_offset(offset, y, dist_1))
  expect_error(check_offset(offset, y, dist_2))
  expect_error(check_offset(offset, y, dist_3), NA)
  expect_error(check_offset(offset, y, dist_4))
  expect_error(check_offset(offset, y, dist_5))
  expect_error(check_offset(offset, y, dist_6))
  expect_error(check_offset(offset, y, dist_7))
  expect_error(check_offset(offset, y, dist_8))
  expect_error(check_offset(offset, y, dist_9))
  expect_error(check_offset(offset, y, dist_10))
  expect_error(check_offset(offset, y, dist_11))
  expect_error(check_offset(offset, y, dist_12))
})


test_that("check_cv_parameters throws an error if cv_folds is not a natural number >= 1", {
  # When cv_folds is not a natural number >=1
  # Then an error is thrown
  expect_error(check_cv_parameters(cv_folds=-1, cv_class_stratify=FALSE, fold_id=NULL,
                                 train_params=training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                                                              shrinkage=0.005, bag_fraction=0.5, id=1:1000, num_train=1000, num_features=6)))
  expect_error(check_cv_parameters(cv_folds=FALSE, cv_class_stratify=FALSE, fold_id=NULL,
                                   train_params=training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                                                                shrinkage=0.005, bag_fraction=0.5, id=1:1000, num_train=1000, num_features=6)))
  expect_error(check_cv_parameters(cv_folds="string", cv_class_stratify=FALSE, fold_id=NULL,
                                   train_params=training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                                                                shrinkage=0.005, bag_fraction=0.5, id=1:1000, num_train=1000, num_features=6)))
  expect_error(check_cv_parameters(cv_folds=c(1, 2), cv_class_stratify=FALSE, fold_id=NULL,
                                   train_params=training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                                                                shrinkage=0.005, bag_fraction=0.5, id=1:1000, num_train=1000, num_features=6)))
  expect_error(check_cv_parameters(cv_folds=2.1, cv_class_stratify=FALSE, fold_id=NULL,
                                   train_params=training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                                                                shrinkage=0.005, bag_fraction=0.5, id=1:1000, num_train=1000, num_features=6)))
  expect_error(check_cv_parameters(cv_folds=NA, cv_class_stratify=FALSE, fold_id=NULL,
                                   train_params=training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                                                                shrinkage=0.005, bag_fraction=0.5, id=1:1000, num_train=1000, num_features=6)))
  expect_error(check_cv_parameters(cv_folds=NaN, cv_class_stratify=FALSE, fold_id=NULL,
                                   train_params=training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                                                                shrinkage=0.005, bag_fraction=0.5, id=1:1000, num_train=1000, num_features=6)))
  expect_error(check_cv_parameters(cv_folds=Inf, cv_class_stratify=FALSE, fold_id=NULL,
                                   train_params=training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                                                                shrinkage=0.005, bag_fraction=0.5, id=1:1000, num_train=1000, num_features=6)))
})

test_that("check_cv_parameters throws an error if cv_class_stratify is not a logical", {
  # When cv_class_stratify is not a logical
  # Then check_cv_parameters throws an error
  expect_error(check_cv_parameters(cv_folds=5, cv_class_stratify="Help", fold_id=NULL,
                                   train_params=training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                                                                shrinkage=0.005, bag_fraction=0.5, id=1:1000, num_train=1000, num_features=6)))
  expect_error(check_cv_parameters(cv_folds=5, cv_class_stratify=NULL, fold_id=NULL,
                                   train_params=training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                                                                shrinkage=0.005, bag_fraction=0.5, id=1:1000, num_train=1000, num_features=6)))
  expect_error(check_cv_parameters(cv_folds=5, cv_class_stratify=NA, fold_id=NULL,
                                   train_params=training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                                                                shrinkage=0.005, bag_fraction=0.5, id=1:1000, num_train=1000, num_features=6)))
  expect_error(check_cv_parameters(cv_folds=5, cv_class_stratify=NaN, fold_id=NULL,
                                   train_params=training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                                                                shrinkage=0.005, bag_fraction=0.5, id=1:1000, num_train=1000, num_features=6)))
  expect_error(check_cv_parameters(cv_folds=5, cv_class_stratify=c(1, 2), fold_id=NULL,
                                   train_params=training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                                                                shrinkage=0.005, bag_fraction=0.5, id=1:1000, num_train=1000, num_features=6)))
  expect_error(check_cv_parameters(cv_folds=5, cv_class_stratify=Inf, fold_id=NULL,
                                   train_params=training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                                                                shrinkage=0.005, bag_fraction=0.5, id=1:1000, num_train=1000, num_features=6)))
})

test_that("check_cv_parameters throws an error if train_params not of class GBMTrainParams", {
  # Given training parameters
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=1:1000, num_train=1000, num_features=6)
  # When class removed
  class(params) <- ""
  
  # Then check_cv_parameters throws an error
  expect_error(check_cv_parameters(cv_folds=5, cv_class_stratify=FALSE, fold_id=NULL,
                                   train_params=params))
})

test_that("check_cv_parameters throws an error if fold_id has observations across different folds", {
  # Given training parameters and obs id
  N <- 2000
  obs_id <- c(rep(1, N/2), rep(2, N/2))
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
         shrinkage=0.005, bag_fraction=0.5, id=obs_id, num_train=N/2, num_features=6)
  
  # When fold_id not NULL and has observations across folds
  cv_folds <- 5
  fold_id <- sample(seq_len(cv_folds), replace=TRUE)
  
  # Then check_cv_parameters throws an error
  expect_error(check_cv_parameters(cv_folds=cv_folds, cv_class_stratify=FALSE, fold_id=fold_id,
                                   train_params=params))
})