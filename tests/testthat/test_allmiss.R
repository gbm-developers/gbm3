context("test all missing vars")

test_that("Columns with all missing stop", {
  
  iris2 <- iris
  iris2$Sepal.Width = NA
  expect_error(gbm(Species == 'setosa' ~ ., data = iris2, distribution = 'bernoulli', n.trees = 2)
               ,'variable\\(s\\) 2\\: Sepal.Width contain only missing values\\.')

  iris2$Petal.Width = NA
  expect_error(gbm(Species == 'setosa' ~ ., data = iris2, distribution = 'bernoulli', n.trees = 2)
               ,'variable\\(s\\) 2, 4\\: Sepal.Width, Petal.Width contain only missing values\\.')
  
  })
