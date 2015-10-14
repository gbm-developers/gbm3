context("Test error messages")

test_that("Columns with all missing stop", {
  
  iris2 <- iris
  iris2$Sepal.Width = NA
  expect_error(
    gbm(Species == 'setosa' ~ ., data = iris2, distribution = 'bernoulli', n.trees = 2)
    ,'variable\\(s\\) 2\\: Sepal.Width contain only missing values\\.'
    )

  iris2$Petal.Width = NA
  expect_error(
    gbm(Species == 'setosa' ~ ., data = iris2, distribution = 'bernoulli', n.trees = 2)
    ,'variable\\(s\\) 2, 4\\: Sepal.Width, Petal.Width contain only missing values\\.'
    )
  
  })


test_that("x longer than y", {
  
  expect_error(
    gbm.fit(x = iris[,c('Sepal.Length', 'Petal.Length')]
            , y = (iris$Species[1:145] == 'setosa')
            , distribution = 'bernoulli', n.trees = 2)
    ,'The number of rows in x does not equal the length of y\\.')

  })


test_that("Excessive levels in x", {
  
  testExcess <- data.frame(
    y = sample(c(0,1), 1025, replace = TRUE)
    ,x1 = runif(1025)
    ,x2 = factor(1:1025)
  )
  
  expect_error(
    gbm.fit(x = testExcess[,c('x1', 'x2')]
            , y = testExcess$y
            , distribution = 'bernoulli', n.trees = 2)
    ,'gbm does not currently handle categorical variables with more than 1024 levels\\. Variable 2\\: x2 has 1025 levels\\.')

  })


test_that("Inacceptable classes", {
  
  testClasses <- data.frame(
    y = sample(c(0,1), 15, replace = TRUE)
    ,x1 = runif(15)
    ,x2 = seq(as.Date('2015-01-01'), as.Date('2015-01-15'), 'days')
  )
  
  expect_error(
    gbm.fit(x = testClasses[,c('x1', 'x2')]
            , y = testClasses$y
            , distribution = 'bernoulli', n.trees = 2)
    ,'variable 2\\: x2 is not of type numeric, ordered, or factor\\.')

  })

