####################
# Author: James Hickey
#
# Series of test to check the functionality that creates strata
#
####################

context("Testing Strata creation:")
test_that("Strata creation function requires GBMData and GBMDist objects", {
  
})

test_that("Distribution object remains unchanged if not CoxPH - strata remain undefined", {
  
})

test_that("Creating strata fills strata, time_order and sorted fields - CoxPH", {
  
})

test_that("Strata can only be created when response is Survival object - CoxPH", {
  
})

test_that("If strata field in distribution object is NULL, all data are put in same strata - CoxPH", {
  
})

test_that("The training responses are sorted according to strata and this order is stored in time_order", {
  
})

test_that("Strata are sorted according to observation id if provided - CoxPH", {
  
})

