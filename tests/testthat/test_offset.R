##############
# Author: James Hickey
# 
# Series of tests on the effect of the offset on the output.
# Help identify if the refactoring is changing the higher level behaviour.
#
##############

test_that("Setting the offset to 0 does not alter the initial value - Adaboost",{
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  
  # random weights if you want to experiment with them
  w <- rexp(N)
  w <- N*w/sum(w)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  
  # Offset 
  offset <- rep(0, N)
  
  # Generate new gbm object
  set.seed(15)
  gbm.no.offset <- gbm(Y~X1+X2+X3,                # formula
                  data=data,                 # dataset
                  weights=w,
                  var.monotone=c(0,0,0),     # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                  distribution="adaboost",
                  n.trees=3000,              # number of trees
                  shrinkage=0.001,           # shrinkage or learning rate, 0.001 to 0.1 usually work
                  interaction.depth=3,       # 1: additive model, 2: two-way interactions, etc
                  bag.fraction = 0.5,        # subsampling fraction, 0.5 is probably best
                  train.fraction = 0.5,      # fraction of data for training, first train.fraction*N used for training
                  cv.folds=5,                # do 5-fold cross-validation
                  n.minobsinnode = 10,       # minimum total weight needed in each node
                  verbose = FALSE)           # don't print progress
  set.seed(15)
  gbm.zero.offset <- gbm(Y~X1+X2+X3,                # formula
                       data=data,                 # dataset
                       weights=w,
                       offset = offset,
                       var.monotone=c(0,0,0),     # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                       distribution="adaboost",
                       n.trees=3000,              # number of trees
                       shrinkage=0.001,           # shrinkage or learning rate, 0.001 to 0.1 usually work
                       interaction.depth=3,       # 1: additive model, 2: two-way interactions, etc
                       bag.fraction = 0.5,        # subsampling fraction, 0.5 is probably best
                       train.fraction = 0.5,      # fraction of data for training, first train.fraction*N used for training
                       cv.folds=5,                # do 5-fold cross-validation
                       n.minobsinnode = 10,       # minimum total weight needed in each node
                       verbose = FALSE)           # don't print progress
  
  expect_equal(gbm.no.offset$initF, gbm.zero.offset$initF)
})

test_that("Setting the offset to 0 does not alter the initial value - Bernoulli",{
  # create some data
  set.seed(1)
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  
  # random weights if you want to experiment with them
  w <- rexp(N)
  w <- N*w/sum(w)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  
  # Offset
  offset <- rep(0, N)
  
  # Generate new gbm object
  set.seed(15)
  gbm.no.offset <- gbm(Y~X1+X2+X3,                # formula
                  data=data,                 # dataset
                  weights=w,
                  var.monotone=c(0,0,0),     # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                  distribution="bernoulli",
                  n.trees=3000,              # number of trees
                  shrinkage=0.001,           # shrinkage or learning rate, 0.001 to 0.1 usually work
                  interaction.depth=3,       # 1: additive model, 2: two-way interactions, etc
                  bag.fraction = 0.5,        # subsampling fraction, 0.5 is probably best
                  train.fraction = 0.5,      # fraction of data for training, first train.fraction*N used for training
                  cv.folds=5,                # do 5-fold cross-validation
                  n.minobsinnode = 10,       # minimum total weight needed in each node
                  verbose = FALSE)           # don't print progress
  set.seed(15)
  gbm.zero.offset <- gbm(Y~X1+X2+X3,                # formula
                       data=data,                 # dataset
                       weights=w,
                       offset=offset,
                       var.monotone=c(0,0,0),     # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                       distribution="bernoulli",
                       n.trees=3000,              # number of trees
                       shrinkage=0.001,           # shrinkage or learning rate, 0.001 to 0.1 usually work
                       interaction.depth=3,       # 1: additive model, 2: two-way interactions, etc
                       bag.fraction = 0.5,        # subsampling fraction, 0.5 is probably best
                       train.fraction = 0.5,      # fraction of data for training, first train.fraction*N used for training
                       cv.folds=5,                # do 5-fold cross-validation
                       n.minobsinnode = 10,       # minimum total weight needed in each node
                       verbose = FALSE)           # don't print progress
  
  expect_equal(gbm.no.offset$initF, gbm.zero.offset$initF)
})

test_that("Setting the offset to 0 does not alter the initial value - CoxPH",{
  require(survival)
  set.seed(1)
  ## GIVEN
  # create some data
  N <- 3000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  f <- 0.5*sin(3*X1 + 5*X2^2 + mu/10)
  tt.surv <- rexp(N,exp(f))
  tt.cens <- rexp(N,0.5)
  delta <- as.numeric(tt.surv <= tt.cens)
  tt <- apply(cbind(tt.surv,tt.cens),1,min)
  
  # throw in some missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  # random weights if you want to experiment with them
  w <- rep(1,N)
  
  data <- data.frame(y=Surv(tt,delta),X1=X1,X2=X2,X3=X3)
  
  #Offset 
  offset <- rep(0, N)
  
  ## WHEN
  # Generate new gbm object
  set.seed(15)
  gbm.no.offset <- gbm(y~X1+X2+X3,       # formula
                  data=data,                 # dataset
                  weights=w,
                  var.monotone=c(0,0,0),     # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                  distribution="coxph",
                  n.trees=3000,              # number of trees
                  shrinkage=0.001,           # shrinkage or learning rate, 0.001 to 0.1 usually work
                  interaction.depth=3,       # 1: additive model, 2: two-way interactions, etc
                  bag.fraction = 0.5,        # subsampling fraction, 0.5 is probably best
                  train.fraction = 0.5,      # fraction of data for training, first train.fraction*N used for training
                  cv.folds = 5,              # do 5-fold cross-validation
                  n.minobsinnode = 10,       # minimum total weight needed in each node
                  keep.data = TRUE,
                  verbose = FALSE)           # don't print progress
  set.seed(15)
  gbm.zero.offset <- gbm(y~X1+X2+X3,       # formula
                       data=data,                 # dataset
                       weights=w,
                       var.monotone=c(0,0,0),     # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                       distribution="coxph",
                       offset=offset,
                       n.trees=3000,              # number of trees
                       shrinkage=0.001,           # shrinkage or learning rate, 0.001 to 0.1 usually work
                       interaction.depth=3,       # 1: additive model, 2: two-way interactions, etc
                       bag.fraction = 0.5,        # subsampling fraction, 0.5 is probably best
                       train.fraction = 0.5,      # fraction of data for training, first train.fraction*N used for training
                       cv.folds = 5,              # do 5-fold cross-validation
                       n.minobsinnode = 10,       # minimum total weight needed in each node
                       keep.data = TRUE,
                       verbose = FALSE)           # don't print progress
  
  expect_equal(gbm.no.offset$initF, gbm.zero.offset$initF)
})

test_that("Setting the offset to 0 does not alter the initial value - Gamma",{
  # Create data
  set.seed(1)
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  X4 <- ordered(sample(letters[1:6],N,replace=T))
  X5 <- factor(sample(letters[1:3],N,replace=T))
  X6 <- 3*runif(N)
  mu <- c(0,0,1,2)[as.numeric(X3)]
  
  SNR <- 10 # signal-to-noise ratio
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  sigma <- sqrt(var(Y)/SNR)
  Y <- Y + abs(rnorm(N,0,sigma))
  
  # create a bunch of missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  w <- rep(1,N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # Offset
  offset <- rep(0, N)
  
  # Generate new gbm object
  set.seed(15)
  gbm.no.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
                  data=data,                   # dataset
                  var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                  distribution="gamma",     # bernoulli, adaboost, gaussian, poisson, coxph, or
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
  set.seed(15)
  gbm.zero.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
                       data=data,                   # dataset
                       var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                       distribution="gamma",     # bernoulli, adaboost, gaussian, poisson, coxph, or
                       # list(name="quantile",alpha=0.05) for quantile regression
                       n.trees=2000,                 # number of trees
                       offset=offset,
                       shrinkage=0.005,             # shrinkage or learning rate, 0.001 to 0.1 usually work
                       interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc
                       bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
                       train.fraction = 0.5,        # fraction of data for training, first train.fraction*N used for training
                       mFeatures = 3,               # Number of features to consider at each node.
                       n.minobsinnode = 10,         # minimum number of obs needed in each node
                       keep.data=TRUE,
                       cv.folds=10,                 # do 10-fold cross-validation
                       verbose = FALSE)             # don't print progress

    expect_equal(gbm.no.offset$initF, gbm.zero.offset$initF)
})

test_that("Setting the offset to 0 does not alter the initial value - Gaussian",{
  # Create data
  set.seed(1)
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
  
  w <- rep(1,N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # Offset
  offset <- rep(0, N)
  
  # Generate new gbm object
  set.seed(15)
  gbm.no.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
                  data=data,                   # dataset
                  var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                  distribution="gaussian",     # bernoulli, adaboost, gaussian, poisson, coxph, or
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
  set.seed(15)
  gbm.zero.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
                       data=data,                   # dataset
                       var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                       distribution="gaussian",     # bernoulli, adaboost, gaussian, poisson, coxph, or
                       # list(name="quantile",alpha=0.05) for quantile regression
                       n.trees=2000,                 # number of trees
                       offset=offset,
                       shrinkage=0.005,             # shrinkage or learning rate, 0.001 to 0.1 usually work
                       interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc
                       bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
                       train.fraction = 0.5,        # fraction of data for training, first train.fraction*N used for training
                       mFeatures = 3,               # Number of features to consider at each node.
                       n.minobsinnode = 10,         # minimum number of obs needed in each node
                       keep.data=TRUE,
                       cv.folds=10,                 # do 10-fold cross-validation
                       verbose = FALSE)             # don't print progress
  
  expect_equal(gbm.no.offset$initF, gbm.zero.offset$initF)
})

test_that("Setting the offset to 0 does not alter the initial value - Laplace",{
  # Create data
  set.seed(1)
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
  
  w <- rep(1,N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # Offset
  offset <- rep(0, N)
  
  # Generate new gbm object
  set.seed(15)
  gbm.no.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
                  data=data,                   # dataset
                  var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                  distribution="laplace",     # bernoulli, adaboost, gaussian, poisson, coxph, or
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
  set.seed(15)
  gbm.zero.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
                       data=data,                   # dataset
                       var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                       distribution="laplace",     # bernoulli, adaboost, gaussian, poisson, coxph, or
                       # list(name="quantile",alpha=0.05) for quantile regression
                       n.trees=2000,                 # number of trees
                       shrinkage=0.005,             # shrinkage or learning rate, 0.001 to 0.1 usually work
                       interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc
                       bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
                       train.fraction = 0.5,        # fraction of data for training, first train.fraction*N used for training
                       mFeatures = 3,               # Number of features to consider at each node.
                       n.minobsinnode = 10,         # minimum number of obs needed in each node
                       keep.data=TRUE,
                       offset=offset,
                       cv.folds=10,                 # do 10-fold cross-validation
                       verbose = FALSE)             # don't print progress
  
  expect_equal(gbm.no.offset$initF, gbm.zero.offset$initF)
})
test_that("Setting the offset to 0 does not alter the initial value - Huberized Hinge Loss",{
  # create some data
  set.seed(1)
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  
  # random weights if you want to experiment with them
  w <- rexp(N)
  w <- N*w/sum(w)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  
  # Offset
  offset <- rep(0, N)
  
  # Generate new gbm object
  set.seed(15)
  gbm.no.offset <- gbm(Y~X1+X2+X3,                # formula
                  data=data,                 # dataset
                  weights=w,
                  var.monotone=c(0,0,0),     # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                  distribution="huberized",
                  n.trees=3000,              # number of trees
                  shrinkage=0.001,           # shrinkage or learning rate, 0.001 to 0.1 usually work
                  interaction.depth=3,       # 1: additive model, 2: two-way interactions, etc
                  bag.fraction = 0.5,        # subsampling fraction, 0.5 is probably best
                  train.fraction = 0.5,      # fraction of data for training, first train.fraction*N used for training
                  cv.folds=5,                # do 5-fold cross-validation
                  n.minobsinnode = 10,       # minimum total weight needed in each node
                  verbose = FALSE)           # don't print progress
  set.seed(15)
  gbm.zero.offset <- gbm(Y~X1+X2+X3,                # formula
                       data=data,                 # dataset
                       weights=w,
                       offset=offset,
                       var.monotone=c(0,0,0),     # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                       distribution="huberized",
                       n.trees=3000,              # number of trees
                       shrinkage=0.001,           # shrinkage or learning rate, 0.001 to 0.1 usually work
                       interaction.depth=3,       # 1: additive model, 2: two-way interactions, etc
                       bag.fraction = 0.5,        # subsampling fraction, 0.5 is probably best
                       train.fraction = 0.5,      # fraction of data for training, first train.fraction*N used for training
                       cv.folds=5,                # do 5-fold cross-validation
                       n.minobsinnode = 10,       # minimum total weight needed in each node
                       verbose = FALSE)           # don't print progress
  
  expect_equal(gbm.no.offset$initF, gbm.zero.offset$initF)
})
test_that("Setting the offset to 0 does not alter the initial value - Pairwise",{
  # create query groups, with an average size of 25 items each
  set.seed(1)
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
  Y.norm <- round(Y/(max(Y) + 0.001))
  
  data.train <- data.frame(Y, query=query, X1, X2, X3)
  data.train.mrr <- data.frame(Y.norm, query=query, X1, X2, X3)
  
  
  # Generate new gbm objects
  set.seed(15)
  gbm.no.offset.ndcg <- gbm(Y~X1+X2+X3,          # formula
                  data=data.train,     # dataset
                  distribution=list(   # loss function:
                    name='pairwise',   # pairwise
                    metric="ndcg",     # ranking metric: normalized discounted cumulative gain
                    group='query'),    # column indicating query groups
                  n.trees=2000,        # number of trees
                  shrinkage=0.005,     # learning rate
                  interaction.depth=3, # number per splits per tree
                  bag.fraction = 0.5,  # subsampling fraction
                  train.fraction = 1,  # fraction of data for training
                  n.minobsinnode = 10, # minimum number of obs for split
                  keep.data=TRUE,      # store copy of input data in model
                  cv.folds=5,          # number of cross validation folds
                  verbose = FALSE,     # don't print progress
                  n.cores = 1)         # use a single core
  
  gbm.no.offset.mrr <- gbm(Y.norm~X1+X2+X3,          # formula
                 data=data.train.mrr,     # dataset
                 distribution=list(   # loss function:
                   name='pairwise',   # pairwise
                   metric="mrr",     # ranking metric: normalized discounted cumulative gain
                   group='query'),    # column indicating query groups
                 n.trees=2000,        # number of trees
                 shrinkage=0.005,     # learning rate
                 interaction.depth=3, # number per splits per tree
                 bag.fraction = 0.5,  # subsampling fraction
                 train.fraction = 1,  # fraction of data for training
                 n.minobsinnode = 10, # minimum number of obs for split
                 keep.data=TRUE,      # store copy of input data in model
                 cv.folds=5,          # number of cross validation folds
                 verbose = FALSE,     # don't print progress
                 n.cores = 1)         # use a single core
  
  gbm.no.offset.map <- gbm(Y.norm~X1+X2+X3,          # formula
                 data=data.train.mrr,     # dataset
                 distribution=list(   # loss function:
                   name='pairwise',   # pairwise
                   metric="map",     # ranking metric: normalized discounted cumulative gain
                   group='query'),    # column indicating query groups
                 n.trees=2000,        # number of trees
                 shrinkage=0.005,     # learning rate
                 interaction.depth=3, # number per splits per tree
                 bag.fraction = 0.5,  # subsampling fraction
                 train.fraction = 1,  # fraction of data for training
                 n.minobsinnode = 10, # minimum number of obs for split
                 keep.data=TRUE,      # store copy of input data in model
                 cv.folds=5,          # number of cross validation folds
                 verbose = FALSE,     # don't print progress
                 n.cores = 1)         # use a single core
  gbm.no.offset.conc <- gbm(Y~X1+X2+X3,          # formula
                  data=data.train,     # dataset
                  distribution=list(   # loss function:
                    name='pairwise',   # pairwise
                    metric="conc",     # ranking metric: concordant pairs
                    group='query'),    # column indicating query groups
                  n.trees=2000,        # number of trees
                  shrinkage=0.005,     # learning rate
                  interaction.depth=3, # number per splits per tree
                  bag.fraction = 0.5,  # subsampling fraction
                  train.fraction = 1,  # fraction of data for training
                  n.minobsinnode = 10, # minimum number of obs for split
                  keep.data=TRUE,      # store copy of input data in model
                  cv.folds=5,          # number of cross validation folds
                  verbose = FALSE,     # don't print progress
                  n.cores = 1)         # use a single core
  set.seed(15)
  gbm.zero.offset.ndcg <- gbm(Y~X1+X2+X3,          # formula
                            data=data.train,     # dataset
                            distribution=list(   # loss function:
                              name='pairwise',   # pairwise
                              metric="ndcg",     # ranking metric: normalized discounted cumulative gain
                              group='query'),    # column indicating query groups
                            n.trees=2000,        # number of trees
                            shrinkage=0.005,     # learning rate
                            interaction.depth=3, # number per splits per tree
                            bag.fraction = 0.5,  # subsampling fraction
                            train.fraction = 1,  # fraction of data for training
                            n.minobsinnode = 10, # minimum number of obs for split
                            keep.data=TRUE,      # store copy of input data in model
                            cv.folds=5,          # number of cross validation folds
                            verbose = FALSE,     # don't print progress
                            offset=rep(0, N),
                            n.cores = 1)         # use a single core
  
  gbm.zero.offset.mrr <- gbm(Y.norm~X1+X2+X3,          # formula
                           data=data.train.mrr,     # dataset
                           distribution=list(   # loss function:
                             name='pairwise',   # pairwise
                             metric="mrr",     # ranking metric: normalized discounted cumulative gain
                             group='query'),    # column indicating query groups
                           n.trees=2000,        # number of trees
                           shrinkage=0.005,     # learning rate
                           interaction.depth=3, # number per splits per tree
                           bag.fraction = 0.5,  # subsampling fraction
                           train.fraction = 1,  # fraction of data for training
                           n.minobsinnode = 10, # minimum number of obs for split
                           keep.data=TRUE,      # store copy of input data in model
                           cv.folds=5,          # number of cross validation folds
                           verbose = FALSE,     # don't print progress
                           offset=rep(0,N),
                           n.cores = 1)         # use a single core
  
  gbm.zero.offset.map <- gbm(Y.norm~X1+X2+X3,          # formula
                           data=data.train.mrr,     # dataset
                           distribution=list(   # loss function:
                             name='pairwise',   # pairwise
                             metric="map",     # ranking metric: normalized discounted cumulative gain
                             group='query'),    # column indicating query groups
                           n.trees=2000,        # number of trees
                           shrinkage=0.005,     # learning rate
                           interaction.depth=3, # number per splits per tree
                           bag.fraction = 0.5,  # subsampling fraction
                           train.fraction = 1,  # fraction of data for training
                           n.minobsinnode = 10, # minimum number of obs for split
                           keep.data=TRUE,      # store copy of input data in model
                           cv.folds=5,          # number of cross validation folds
                           verbose = FALSE,     # don't print progress
                           offset=rep(0,N),
                           n.cores = 1)         # use a single core
  gbm.zero.offset.conc <- gbm(Y~X1+X2+X3,          # formula
                            data=data.train,     # dataset
                            distribution=list(   # loss function:
                              name='pairwise',   # pairwise
                              metric="conc",     # ranking metric: concordant pairs
                              group='query'),    # column indicating query groups
                            n.trees=2000,        # number of trees
                            shrinkage=0.005,     # learning rate
                            interaction.depth=3, # number per splits per tree
                            bag.fraction = 0.5,  # subsampling fraction
                            train.fraction = 1,  # fraction of data for training
                            n.minobsinnode = 10, # minimum number of obs for split
                            keep.data=TRUE,      # store copy of input data in model
                            cv.folds=5,          # number of cross validation folds
                            verbose = FALSE,     # don't print progress
                            offset = rep(0,N),
                            n.cores = 1)         # use a single core
  
  expect_equal(gbm.no.offset.conc$initF, gbm.zero.offset.conc$initF)
  expect_equal(gbm.no.offset.mrr$initF, gbm.zero.offset.mrr$initF)
  expect_equal(gbm.no.offset.map$initF, gbm.zero.offset.map$initF)
  expect_equal(gbm.no.offset.ndcg$initF, gbm.zero.offset.ndcg$initF)
})
test_that("Setting the offset to 0 does not alter the initial value - Poisson",{
  # Create data
  set.seed(1)
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  X4 <- ordered(sample(letters[1:6],N,replace=T))
  X5 <- factor(sample(letters[1:3],N,replace=T))
  X6 <- 3*runif(N)
  mu <- c(0,0,1,2)[as.numeric(X3)]
  
  SNR <- 10 # signal-to-noise ratio
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  sigma <- sqrt(var(Y)/SNR)
  Y <- round(Y + abs(rnorm(N,0,sigma))) # Ensure it is normal
  
  # create a bunch of missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  w <- rep(1,N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # Offset
  offset <- rep(0, N)
  
  # Generate new gbm object
  set.seed(15)
  gbm.no.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
                  data=data,                   # dataset
                  var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                  distribution="poisson",     # bernoulli, adaboost, gaussian, poisson, coxph, or
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
  set.seed(15)
  gbm.zero.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
                       data=data,                   # dataset
                       var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                       distribution="poisson",     # bernoulli, adaboost, gaussian, poisson, coxph, or
                       # list(name="quantile",alpha=0.05) for quantile regression
                       n.trees=2000,                 # number of trees
                       offset=offset,
                       shrinkage=0.005,             # shrinkage or learning rate, 0.001 to 0.1 usually work
                       interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc
                       bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
                       train.fraction = 0.5,        # fraction of data for training, first train.fraction*N used for training
                       mFeatures = 3,               # Number of features to consider at each node.
                       n.minobsinnode = 10,         # minimum number of obs needed in each node
                       keep.data=TRUE,
                       cv.folds=10,                 # do 10-fold cross-validation
                       verbose = FALSE)             # don't print progress
  
  expect_equal(gbm.no.offset$initF, gbm.zero.offset$initF)
})
test_that("Setting the offset to 0 does not alter the initial value - Quantile Reg",{
  # Create data
  set.seed(1)
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
  
  w <- rep(1,N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  #Offset
  offset <- rep(0, N)
  
  # Generate new gbm object
  set.seed(15)
  gbm.no.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
                  data=data,                   # dataset
                  var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                  distribution= list(name="quantile",alpha=0.95), # for quantile regression
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
  set.seed(15)
  gbm.zero.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
                       data=data,                   # dataset
                       var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                       distribution= list(name="quantile",alpha=0.95), # for quantile regression
                       n.trees=2000,                 # number of trees
                       shrinkage=0.005,             # shrinkage or learning rate, 0.001 to 0.1 usually work
                       interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc
                       bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
                       train.fraction = 0.5,        # fraction of data for training, first train.fraction*N used for training
                       mFeatures = 3,               # Number of features to consider at each node.
                       n.minobsinnode = 10,         # minimum number of obs needed in each node
                       keep.data=TRUE,
                       offset=offset,
                       cv.folds=10,                 # do 10-fold cross-validation
                       verbose = FALSE)             # don't print progress
  
  expect_equal(gbm.no.offset$initF, gbm.zero.offset$initF)
})
test_that("Setting the offset to 0 does not alter the initial value - T dist",{
  # Create data
  set.seed(1)
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
  
  w <- rep(1,N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # Offset
  offset <- rep(0, N)
  
  # Generate new gbm object
  set.seed(15)
  gbm.no.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
                  data=data,                   # dataset
                  var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                  distribution="tdist",     # bernoulli, adaboost, gaussian, poisson, coxph, or
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
  set.seed(15)
  gbm.zero.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
                       data=data,                   # dataset
                       var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                       distribution="tdist",     # bernoulli, adaboost, gaussian, poisson, coxph, or
                       # list(name="quantile",alpha=0.05) for quantile regression
                       n.trees=2000,                 # number of trees
                       shrinkage=0.005,             # shrinkage or learning rate, 0.001 to 0.1 usually work
                       interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc
                       bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
                       train.fraction = 0.5,        # fraction of data for training, first train.fraction*N used for training
                       mFeatures = 3,               # Number of features to consider at each node.
                       n.minobsinnode = 10,         # minimum number of obs needed in each node
                       keep.data=TRUE,
                       offset=offset,
                       cv.folds=10,                 # do 10-fold cross-validation
                       verbose = FALSE)             # don't print progress
  
  expect_equal(gbm.no.offset$initF, gbm.zero.offset$initF)
})
test_that("Setting the offset to 0 does not alter the initial value - Tweedie",{
  # Create data
  set.seed(1)
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  X4 <- ordered(sample(letters[1:6],N,replace=T))
  X5 <- factor(sample(letters[1:3],N,replace=T))
  X6 <- 3*runif(N)
  mu <- c(0,0,1,2)[as.numeric(X3)]
  
  SNR <- 10 # signal-to-noise ratio
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  sigma <- sqrt(var(Y)/SNR)
  Y <- Y + abs(rnorm(N,0,sigma))
  
  # create a bunch of missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  w <- rep(1,N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # Offset
  offset <- rep(0, N)
  
  # Generate new gbm object
  set.seed(15)
  gbm.no.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
                  data=data,                   # dataset
                  var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                  distribution="tweedie",     # bernoulli, adaboost, gaussian, poisson, coxph, or
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
  set.seed(15)
  gbm.zero.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
                       data=data,                   # dataset
                       var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                       distribution="tweedie",     # bernoulli, adaboost, gaussian, poisson, coxph, or
                       # list(name="quantile",alpha=0.05) for quantile regression
                       n.trees=2000,                 # number of trees
                       shrinkage=0.005,             # shrinkage or learning rate, 0.001 to 0.1 usually work
                       interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc
                       bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
                       train.fraction = 0.5,        # fraction of data for training, first train.fraction*N used for training
                       mFeatures = 3,               # Number of features to consider at each node.
                       n.minobsinnode = 10,         # minimum number of obs needed in each node
                       keep.data=TRUE,
                       offset=offset,
                       cv.folds=10,                 # do 10-fold cross-validation
                       verbose = FALSE)             # don't print progress
  
  expect_equal(gbm.no.offset$initF, gbm.zero.offset$initF)
})
test_that("Increasing the offset reduces the initial value - Adaboost", {
  # create some data
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
  Y <- rbinom(N,1,p)
  
  # random weights if you want to experiment with them
  w <- rexp(N)
  w <- N*w/sum(w)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
  
  # Offset 
  offset.smaller <- rep(1, N)
  offset.larger <- rep(20, N)
  
  # Generate new gbm object
  set.seed(15)
  gbm.smaller.offset <- gbm(Y~X1+X2+X3,                # formula
                       data=data,                 # dataset
                       weights=w,
                       offset=offset.smaller,
                       var.monotone=c(0,0,0),     # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                       distribution="adaboost",
                       n.trees=3000,              # number of trees
                       shrinkage=0.001,           # shrinkage or learning rate, 0.001 to 0.1 usually work
                       interaction.depth=3,       # 1: additive model, 2: two-way interactions, etc
                       bag.fraction = 0.5,        # subsampling fraction, 0.5 is probably best
                       train.fraction = 0.5,      # fraction of data for training, first train.fraction*N used for training
                       cv.folds=5,                # do 5-fold cross-validation
                       n.minobsinnode = 10,       # minimum total weight needed in each node
                       verbose = FALSE)           # don't print progress
  set.seed(15)
  gbm.larger.offset <- gbm(Y~X1+X2+X3,                # formula
                         data=data,                 # dataset
                         weights=w,
                         offset = offset.larger,
                         var.monotone=c(0,0,0),     # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                         distribution="adaboost",
                         n.trees=3000,              # number of trees
                         shrinkage=0.001,           # shrinkage or learning rate, 0.001 to 0.1 usually work
                         interaction.depth=3,       # 1: additive model, 2: two-way interactions, etc
                         bag.fraction = 0.5,        # subsampling fraction, 0.5 is probably best
                         train.fraction = 0.5,      # fraction of data for training, first train.fraction*N used for training
                         cv.folds=5,                # do 5-fold cross-validation
                         n.minobsinnode = 10,       # minimum total weight needed in each node
                         verbose = FALSE)           # don't print progress
  
  expect_true(gbm.smaller.offset$initF - gbm.larger.offset$initF > 0)
})
test_that("Increasing the offset reduces the initial value - Poisson", {
  # Create data
  set.seed(1)
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  X4 <- ordered(sample(letters[1:6],N,replace=T))
  X5 <- factor(sample(letters[1:3],N,replace=T))
  X6 <- 3*runif(N)
  mu <- c(0,0,1,2)[as.numeric(X3)]
  
  SNR <- 10 # signal-to-noise ratio
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  sigma <- sqrt(var(Y)/SNR)
  Y <- round(Y + abs(rnorm(N,0,sigma))) # Ensure it is normal
  
  # create a bunch of missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  w <- rep(1,N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # Offset 
  offset.smaller <- rep(1, N)
  offset.larger <- rep(20, N)
  
  # Generate new gbm object
  set.seed(15)
  gbm.smaller.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
                            data=data,                   # dataset
                            var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                            distribution="poisson",     # bernoulli, adaboost, gaussian, poisson, coxph, or
                            # list(name="quantile",alpha=0.05) for quantile regression
                            n.trees=2000,                 # number of trees
                            offset=offset.smaller,
                            shrinkage=0.005,             # shrinkage or learning rate, 0.001 to 0.1 usually work
                            interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc
                            bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
                            train.fraction = 0.5,        # fraction of data for training, first train.fraction*N used for training
                            mFeatures = 3,               # Number of features to consider at each node.
                            n.minobsinnode = 10,         # minimum number of obs needed in each node
                            keep.data=TRUE,
                            cv.folds=10,                 # do 10-fold cross-validation
                            verbose = FALSE)             # don't print progress
  set.seed(15)
  gbm.larger.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
                         data=data,                   # dataset
                         var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                         distribution="poisson",     # bernoulli, adaboost, gaussian, poisson, coxph, or
                         # list(name="quantile",alpha=0.05) for quantile regression
                         n.trees=2000,                 # number of trees
                         offset=offset.larger,
                         shrinkage=0.005,             # shrinkage or learning rate, 0.001 to 0.1 usually work
                         interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc
                         bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
                         train.fraction = 0.5,        # fraction of data for training, first train.fraction*N used for training
                         mFeatures = 3,               # Number of features to consider at each node.
                         n.minobsinnode = 10,         # minimum number of obs needed in each node
                         keep.data=TRUE,
                         cv.folds=10,                 # do 10-fold cross-validation
                         verbose = FALSE)             # don't print progress
    
  
  expect_true(gbm.smaller.offset$initF - gbm.larger.offset$initF > 0)
})