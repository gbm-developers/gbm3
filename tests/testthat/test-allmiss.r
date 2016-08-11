context("Test error messages - new api")
test_that("Columns with all missing stop", {
  
  iris2 <- iris
  iris2$Sepal.Width = NA
  expect_error(
    gbmt(Species == 'setosa' ~ ., data = iris2, distribution = gbm_dist('Bernoulli'), train_params = training_params(num_trees=2))
    ,'variable 2: Sepal.Width is not of type - numeric, ordered or factor.'
  )
  
  iris2$Petal.Width = NA
  expect_error(
    gbmt(Species == 'setosa' ~ ., data = iris2, distribution = gbm_dist('Bernoulli'), train_params = training_params(num_trees=2))
    ,'variable 2, 4: Sepal.Width, Petal.Width is not of type - numeric, ordered or factor.'
  )
  
})

test_that("x longer than y", {
  
  expect_error(
    gbmt_fit(x = iris[,c('Sepal.Length', 'Petal.Length')]
            , y = (iris$Species[1:145] == 'setosa')
            , distribution = gbm_dist('Bernoulli'), train_params=training_params(num_trees = 2))
    ,'The number of rows in x does not equal the length of y\\.')
  
})

test_that("Excessive levels in x", {
  
  testExcess <- data.frame(
    y = sample(c(0,1), 1025, replace = TRUE)
    ,x1 = runif(1025)
    ,x2 = factor(1:1025)
  )
  
  expect_error(
    gbmt_fit(x = testExcess[,c('x1', 'x2')]
            , y = testExcess$y
            , distribution = gbm_dist("Bernoulli"), train_params=training_params(num_trees = 2))
    ,'gbm does not currently handle categorical variables with more than 1024 levels\\. Variable 2\\: x2 has 1025 levels\\.')
  
})

test_that("Inacceptable classes", {
  
  testClasses <- data.frame(
    y = sample(c(0,1), 15, replace = TRUE)
    ,x1 = runif(15)
    ,x2 = seq(as.Date('2015-01-01'), as.Date('2015-01-15'), 'days')
  )
  
  expect_error(
    gbmt_fit(x = testClasses[,c('x1', 'x2')]
            , y = testClasses$y
            , distribution = gbm_dist('Bernoulli'), train_params=training_params(num_trees = 2, min_num_obs_in_node = 1))
    ,'variable 2\\: x2 is not of type - numeric, ordered or factor\\.')
  
})

context("Test error messages - old api")
test_that("Columns with all missing stop", {
  
  iris2 <- iris
  iris2$Sepal.Width = NA
  expect_error(
    gbm(Species == 'setosa' ~ ., data = iris2, distribution = 'Bernoulli', n.trees = 2)
    ,'variable 2: Sepal.Width is not of type - numeric, ordered or factor.'
    )

  iris2$Petal.Width = NA
  expect_error(
    gbm(Species == 'setosa' ~ ., data = iris2, distribution = 'Bernoulli', n.trees = 2)
    ,'variable 2, 4: Sepal.Width, Petal.Width is not of type - numeric, ordered or factor.'
    )
  
})

test_that("x longer than y", {
  
  expect_error(
    gbm.fit(x = iris[,c('Sepal.Length', 'Petal.Length')]
            , y = (iris$Species[1:145] == 'setosa')
            , distribution = 'Bernoulli', n.trees = 2)
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
            , distribution = 'Bernoulli', n.trees = 2)
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
            , distribution = 'Bernoulli', n.trees = 2, n.minobsinnode = 1)
    ,'variable 2\\: x2 is not of type - numeric, ordered or factor\\.')

})

