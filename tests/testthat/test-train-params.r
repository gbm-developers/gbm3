####################
# Author: James Hickey
#
# Series of test to validate the GBMTrainParams objects.
#
####################

##### Definition #####
context("Check the definition of the GBMTrainParams object")
test_that("Training parameters are of class GBMTrainParams", {
  expect_true("GBMTrainParams" %in% class(training_params(num_train=100, id=seq_len(100), num_features = 3)))
})

test_that("GBMTrainParams structure is correct", {
  expect_equal(names(training_params(num_train=100, id=seq_len(100), num_features=3)), 
               c("num_trees", "interaction_depth", "min_num_obs_in_node", "shrinkage",
                 "bag_fraction", "id", "num_train", "num_features", "num_rows_per_obs"))
})

test_that("Calculation of number of rows per observation is correct", {
  # Given training obs and ids
  train_no1 <- 100
  train_no2 <- 25
  id1 <- seq_len(100)
  id2 <- rep(1:25, 4)
  
  # Calculates number of rows per observation correctly on construction
  expect_equal(training_params(num_train = train_no1, id=id1, num_features = 3)$num_rows_per_obs, table(id=id1[order(id1)]))
  expect_equal(training_params(num_train = train_no2, id=id2, num_features = 3)$num_rows_per_obs, table(id=id2[order(id2)]))
})


##### Errors #####
context("Using invalid contruction parameters will thrown an error")
test_that("Number of trees must be a positive integer", {
  expect_error(training_params(num_trees=-0.1, num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(num_trees=-1, num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(num_trees=c(1, 2), num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(num_trees=Inf, num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(num_trees=NULL, num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(num_trees=NA, num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(num_trees="Incorrect", num_train=100, id=seq_len(100), num_features=3))
})

test_that("Interaction depth must be a positive integer", {
  expect_error(training_params(interaction_depth=-0.1, num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(interaction_depth=-1, num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(interaction_depth=c(1, 2), num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(interaction_depth=Inf, num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(interaction_depth=NULL, num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(interaction_depth=NA, num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(interaction_depth="Incorrect", num_train=100, id=seq_len(100), num_features=3))
})

test_that("Minimum number of node observations must be a positive integer", {
  expect_error(training_params(min_num_obs_in_node=-0.1, num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(min_num_obs_in_node=-1, num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(min_num_obs_in_node=c(1, 2), num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(min_num_obs_in_node=Inf, num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(min_num_obs_in_node=NULL, num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(min_num_obs_in_node=NA, num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(min_num_obs_in_node="Incorrect", num_train=100, id=seq_len(100), num_features=3))
})

test_that("Number of training rows must be a positive integer", {
  expect_error(training_params(num_train=-0.1, id=seq_len(100), num_features=3))
  expect_error(training_params(num_train=-1,  id=seq_len(100), num_features=3))
  expect_error(training_params(num_train=c(1, 2), id=seq_len(100), num_features=3))
  expect_error(training_params(num_train=Inf,  id=seq_len(100), num_features=3))
  expect_error(training_params(num_train=NULL,  id=seq_len(100), num_features=3))
  expect_error(training_params(num_train=NA,  id=seq_len(100), num_features=3))
  expect_error(training_params(num_train="Incorrect",  id=seq_len(100), num_features=3))
})

test_that("Number of features for tree growing must be a positive integer", {
  expect_error(training_params(num_train=100, id=seq_len(100), num_features=-0.1))
  expect_error(training_params(num_train=100,  id=seq_len(100), num_features=-1))
  expect_error(training_params(num_train=100, id=seq_len(100), num_features=c(1,2)))
  expect_error(training_params(num_train=100,  id=seq_len(100), num_features=Inf))
  expect_error(training_params(num_train=100,  id=seq_len(100), num_features=NULL))
  expect_error(training_params(num_train=100,  id=seq_len(100), num_features=NA))
  expect_error(training_params(num_train=100,  id=seq_len(100), num_features="Incorrect"))
})

test_that("Shrinkage must be a double", {
  expect_error(training_params(shrinkage=c(1, 2), num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(shrinkage=Inf, num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(shrinkage=NULL, num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(shrinkage=NA, num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(shrinkage="Incorrect", num_train=100, id=seq_len(100), num_features=3))
})

test_that("Bag fraction must be a double between 0.0 and 1.0", {
  expect_error(training_params(bag_fraction=-0.1, num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(bag_fraction=-1, num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(bag_fraction=1.5, num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(bag_fraction=c(1, 2), num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(bag_fraction=Inf, num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(bag_fraction=NULL, num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(bag_fraction=NA, num_train=100, id=seq_len(100), num_features=3))
  expect_error(training_params(bag_fraction="Incorrect", num_train=100, id=seq_len(100), num_features=3))
})

test_that("Id must be a vector of integers", {
  expect_error(training_params(num_train=100, id=-0.1, num_features=3))
  expect_error(training_params(num_train=100, id=c(0.1, 0.2, 3), num_features=3))
  expect_error(training_params(num_train=100, id=Inf, num_features=3))
  expect_error(training_params(num_train=100, id="Incorrect", num_features=3))
  expect_error(training_params(num_train=100, id=NA, num_features=3))
  expect_error(training_params(num_train=100, id=NULL, num_features=3))
})

test_that("Too few observations throws error", {
  expect_error(training_params(num_train=1, id=1, num_features=3))
})

test_that("Number of training observations must not exceed maximum amount id'ed", {
  expect_error(training_params(num_train=100, id=rep(1:25, 4), num_feactures=3))
})
