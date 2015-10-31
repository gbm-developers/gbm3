context("test ir measures")

test_that("area under curve", {
  
  AUC1 <- gbm.roc.area(c(0,1,0,1,0,1), c(.5,.85,.5,.85,.25,.25))
  expect_true(round(AUC1, 4) == .7222)
  
  AUC2 <- gbm.roc.area(c(rep(1, 10), rep(0, 10))
    , c(0.66, 0.63, 0.44, 0.83, 0.22, 0.33, 0.98, 0.45, 0.23, 0.61
        ,0.73, 0.03, 0.53, 0.58, 0.46, 0.47, 0.22, 0.45, 0.12, 0.42))
  expect_true(AUC2 == .64)
  
})
