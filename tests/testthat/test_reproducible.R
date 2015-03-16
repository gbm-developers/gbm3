context("reproducibility")

test_that("Setting the seed causes result to be reproducible", {
  set.seed(18900217)
  mod <- gbm(Species ~ ., data=iris, distribution="multinomial", cv.folds=3, n.trees=100, shrinkage=.1)
  nt1 <- gbm.perf(mod, method="cv", plot.it=FALSE)
  ri1 <- relative.influence(mod, n.trees=nt1)

  set.seed(18900217)
  mod <- gbm(Species ~ ., data=iris, distribution="multinomial", cv.folds=3, n.trees=100, shrinkage=.1)
  nt2 <- gbm.perf(mod, method="cv", plot.it=FALSE)
  ri2 <- relative.influence(mod, n.trees=nt2)
  
  expect_that(nt1, equals(nt2), label="Number of trees match when same seed is used")
  expect_that(ri1, equals(ri2), label="Relative influences match when same seed is used")
})

test_that("Setting different seeds causes result to be different", {
  set.seed(18900217)
  mod <- gbm(Species ~ ., data=iris, distribution="multinomial", cv.folds=3, n.trees=100, shrinkage=.1)
  nt1 <- gbm.perf(mod, method="cv", plot.it=FALSE)
  ri1 <- relative.influence(mod, n.trees=nt1)

  set.seed(19620729)
  mod <- gbm(Species ~ ., data=iris, distribution="multinomial", cv.folds=3, n.trees=100, shrinkage=.1)
  nt2 <- gbm.perf(mod, method="cv", plot.it=FALSE)
  ri2 <- relative.influence(mod, n.trees=nt2)

  expect_that(nt1, not(equals(nt2)), label="Number of trees don't match when different seeds are used")
  expect_that(ri1, not(equals(ri2)), label="Relative influences don't match when different seeds are used")
})

