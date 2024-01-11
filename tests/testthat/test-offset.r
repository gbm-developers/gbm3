##############
# Author: James Hickey
# 
# Series of tests on the effect of the offset on the output.
# Help identify if the refactoring is changing the higher level behaviour.
#
##############
context("Testing 0 offset - old API")
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
  gbm.no.offset <- gbm(Y~X1+X2+X3,
                  data=data,
                  weights=w,
                  var.monotone=c(0,0,0),
                  distribution="adaboost",
                  n.trees=100,
                  shrinkage=0.001,
                  interaction.depth=3,
                  bag.fraction = 0.5,
                  train.fraction = 0.5,
                  cv.folds=1,
                  n.minobsinnode = 10,
                  verbose = FALSE)
  set.seed(15)
  gbm.zero.offset <- gbm(Y~X1+X2+X3,
                       data=data,
                       weights=w,
                       offset = offset,
                       var.monotone=c(0,0,0),
                       distribution="adaboost",
                       n.trees=100,
                       shrinkage=0.001,
                       interaction.depth=3,
                       bag.fraction = 0.5,
                       train.fraction = 0.5,
                       cv.folds=1,
                       n.minobsinnode = 10,
                       verbose = FALSE)
  
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
  gbm.no.offset <- gbm(Y~X1+X2+X3,
                  data=data,
                  weights=w,
                  var.monotone=c(0,0,0),
                  distribution="bernoulli",
                  n.trees=100,
                  shrinkage=0.001,
                  interaction.depth=3,
                  bag.fraction = 0.5,
                  train.fraction = 0.5,
                  cv.folds=1,
                  n.minobsinnode = 10,
                  verbose = FALSE)
  set.seed(15)
  gbm.zero.offset <- gbm(Y~X1+X2+X3,
                       data=data,
                       weights=w,
                       offset=offset,
                       var.monotone=c(0,0,0),
                       distribution="bernoulli",
                       n.trees=100,
                       shrinkage=0.001,
                       interaction.depth=3,
                       bag.fraction = 0.5,
                       train.fraction = 0.5,
                       cv.folds=1,
                       n.minobsinnode = 10,
                       verbose = FALSE)
  
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
  gbm.no.offset <- gbm(y~X1+X2+X3,
                  data=data,
                  weights=w,
                  var.monotone=c(0,0,0),
                  distribution="coxph",
                  n.trees=100,
                  shrinkage=0.001,
                  interaction.depth=3,
                  bag.fraction = 0.5,
                  train.fraction = 0.5,
                  cv.folds = 1,
                  n.minobsinnode = 10,
                  keep.data = TRUE,
                  verbose = FALSE)
  set.seed(15)
  gbm.zero.offset <- gbm(y~X1+X2+X3,
                       data=data,
                       weights=w,
                       var.monotone=c(0,0,0),
                       distribution="CoxPH",
                       offset=offset,
                       n.trees=100,
                       shrinkage=0.001,
                       interaction.depth=3,
                       bag.fraction = 0.5,
                       train.fraction = 0.5,
                       cv.folds = 1,
                       n.minobsinnode = 10,
                       keep.data = TRUE,
                       verbose = FALSE)
  
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
  gbm.no.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,
                  data=data,
                  var.monotone=c(0,0,0,0,0,0),
                  distribution="gamma",
                  n.trees=100,
                  shrinkage=0.005,
                  interaction.depth=3,
                  bag.fraction = 0.5,
                  train.fraction = 0.5,
                  mFeatures = 3,
                  n.minobsinnode = 10,
                  keep.data=TRUE,
                  cv.folds=1,
                  verbose = FALSE)
  set.seed(15)
  gbm.zero.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,
                       data=data,
                       var.monotone=c(0,0,0,0,0,0),
                       distribution="Gamma",
                       n.trees=100,
                       offset=offset,
                       shrinkage=0.005,
                       interaction.depth=3,
                       bag.fraction = 0.5,
                       train.fraction = 0.5,
                       mFeatures = 3,
                       n.minobsinnode = 10,
                       keep.data=TRUE,
                       cv.folds=1,
                       verbose = FALSE)

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
  gbm.no.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,
                  data=data,
                  var.monotone=c(0,0,0,0,0,0),
                  distribution="gaussian",
                  n.trees=100,
                  shrinkage=0.005,
                  interaction.depth=3,
                  bag.fraction = 0.5,
                  train.fraction = 0.5,
                  mFeatures = 3,
                  n.minobsinnode = 10,
                  keep.data=TRUE,
                  cv.folds=1,
                  verbose = FALSE)
  set.seed(15)
  gbm.zero.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,
                       data=data,
                       var.monotone=c(0,0,0,0,0,0),
                       distribution="gaussian",
                       n.trees=100,
                       offset=offset,
                       shrinkage=0.005,
                       interaction.depth=3,
                       bag.fraction = 0.5,
                       train.fraction = 0.5,
                       mFeatures = 3,
                       n.minobsinnode = 10,
                       keep.data=TRUE,
                       cv.folds=1,
                       verbose = FALSE)
  
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
  gbm.no.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,
                  data=data,
                  var.monotone=c(0,0,0,0,0,0),
                  distribution="laplace",
                  n.trees=100,
                  shrinkage=0.005,
                  interaction.depth=3,
                  bag.fraction = 0.5,
                  train.fraction = 0.5,
                  mFeatures = 3,
                  n.minobsinnode = 10,
                  keep.data=TRUE,
                  cv.folds=1,
                  verbose = FALSE)
  set.seed(15)
  gbm.zero.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,
                       data=data,
                       var.monotone=c(0,0,0,0,0,0),
                       distribution="laplace",
                       n.trees=100,
                       shrinkage=0.005,
                       interaction.depth=3,
                       bag.fraction = 0.5,
                       train.fraction = 0.5,
                       mFeatures = 3,
                       n.minobsinnode = 10,
                       keep.data=TRUE,
                       offset=offset,
                       cv.folds=1,
                       verbose = FALSE)
  
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
  gbm.no.offset <- gbm(Y~X1+X2+X3,
                  data=data,
                  weights=w,
                  var.monotone=c(0,0,0),
                  distribution="huberized",
                  n.trees=100,
                  shrinkage=0.001,
                  interaction.depth=3,
                  bag.fraction = 0.5,
                  train.fraction = 0.5,
                  cv.folds=1,
                  n.minobsinnode = 10,
                  verbose = FALSE)
  set.seed(15)
  gbm.zero.offset <- gbm(Y~X1+X2+X3,
                       data=data,
                       weights=w,
                       offset=offset,
                       var.monotone=c(0,0,0),
                       distribution="huberized",
                       n.trees=100,
                       shrinkage=0.001,
                       interaction.depth=3,
                       bag.fraction = 0.5,
                       train.fraction = 0.5,
                       cv.folds=1,
                       n.minobsinnode = 10,
                       verbose = FALSE)
  
  expect_equal(gbm.no.offset$initF, gbm.zero.offset$initF)
})

test_that("Setting the offset to 0 does not alter the initial value - Pairwise",{
  skip("Skipping pairwise")
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
  gbm.no.offset.ndcg <- gbm(Y~X1+X2+X3,
                  data=data.train,     
                  distribution=list(   
                    name='pairwise',   
                    metric="ndcg",     
                    group="query"),    
                  n.trees=100,        
                  shrinkage=0.005,     
                  interaction.depth=3, 
                  bag.fraction = 0.5,  
                  train.fraction = 1,  
                  n.minobsinnode = 10, 
                  keep.data=TRUE,      
                  cv.folds=1,          
                  verbose = FALSE     
                  )         
  
  gbm.no.offset.mrr <- gbm(Y.norm~X1+X2+X3, 
                 data=data.train.mrr,     
                 distribution=list(  
                   name='pairwise',  
                   metric="mrr",
                   group='query'),   
                 n.trees=100,        
                 shrinkage=0.005,    
                 interaction.depth=3,
                 bag.fraction = 0.5, 
                 train.fraction = 1, 
                 n.minobsinnode = 10,
                 keep.data=TRUE,     
                 cv.folds=1,         
                 verbose = FALSE     
                 )         
  
  gbm.no.offset.map <- gbm(Y.norm~X1+X2+X3,
                 data=data.train.mrr,
                 distribution=list(  
                   name='pairwise',  
                   metric="map",     
                   group='query'),   
                 n.trees=100,        
                 shrinkage=0.005,    
                 interaction.depth=3,
                 bag.fraction = 0.5, 
                 train.fraction = 1, 
                 n.minobsinnode = 10,
                 keep.data=TRUE,     
                 cv.folds=1,         
                 verbose = FALSE     
                 )         
  gbm.no.offset.conc <- gbm(Y~X1+X2+X3,  
                  data=data.train,    
                  distribution=list(  
                    name='pairwise',  
                    metric="conc",    
                    group='query'),   
                  n.trees=100,        
                  shrinkage=0.005,    
                  interaction.depth=3,
                  bag.fraction = 0.5, 
                  train.fraction = 1, 
                  n.minobsinnode = 10,
                  keep.data=TRUE,     
                  cv.folds=1,         
                  verbose = FALSE     
                  )
  set.seed(15)
  gbm.zero.offset.ndcg <- gbm(Y~X1+X2+X3,       
                            data=data.train,    
                            distribution=list(  
                              name='pairwise',  
                              metric="ndcg",     
                              group='query'),   
                            n.trees=100,        
                            shrinkage=0.005,    
                            interaction.depth=3,
                            bag.fraction = 0.5, 
                            train.fraction = 1, 
                            n.minobsinnode = 10,
                            keep.data=TRUE,     
                            cv.folds=1,         
                            verbose = FALSE,    
                            offset=rep(0, N)
                            )
  
  gbm.zero.offset.mrr <- gbm(Y.norm~X1+X2+X3,  
                           data=data.train.mrr,
                           distribution=list(  
                             name='pairwise',  
                             metric="mrr",     
                             group='query'),   
                           n.trees=100,        
                           shrinkage=0.005,    
                           interaction.depth=3,
                           bag.fraction = 0.5, 
                           train.fraction = 1, 
                           n.minobsinnode = 10,
                           keep.data=TRUE,     
                           cv.folds=1,         
                           verbose = FALSE,    
                           offset=rep(0,N)
                           )         
  
  gbm.zero.offset.map <- gbm(Y.norm~X1+X2+X3,  
                           data=data.train.mrr,
                           distribution=list(  
                             name='pairwise',  
                             metric="map",     
                             group='query'),   
                           n.trees=100,        
                           shrinkage=0.005,    
                           interaction.depth=3,
                           bag.fraction = 0.5, 
                           train.fraction = 1, 
                           n.minobsinnode = 10,
                           keep.data=TRUE,     
                           cv.folds=1,         
                           verbose = FALSE,    
                           offset=rep(0,N)
                           )         
  gbm.zero.offset.conc <- gbm(Y~X1+X2+X3,       
                            data=data.train,    
                            distribution=list(  
                              name='pairwise',  
                              metric="conc",    
                              group='query'),   
                            n.trees=100,        
                            shrinkage=0.005,    
                            interaction.depth=3,
                            bag.fraction = 0.5, 
                            train.fraction = 1, 
                            n.minobsinnode = 10,
                            keep.data=TRUE,     
                            cv.folds=1,         
                            verbose = FALSE,    
                            offset = rep(0,N))
  
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
  gbm.no.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,   
                  data=data,                  
                  var.monotone=c(0,0,0,0,0,0),
                  distribution="poisson",     
                  n.trees=100,                
                  shrinkage=0.005,            
                  interaction.depth=3,        
                  bag.fraction = 0.5,         
                  train.fraction = 0.5,       
                  mFeatures = 3,              
                  n.minobsinnode = 10,        
                  keep.data=TRUE,
                  cv.folds=1,                 
                  verbose = FALSE)            
  set.seed(15)
  gbm.zero.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,      
                       data=data,                  
                       var.monotone=c(0,0,0,0,0,0),
                       distribution="poisson",     
                       n.trees=100,                
                       offset=offset,
                       shrinkage=0.005,            
                       interaction.depth=3,        
                       bag.fraction = 0.5,         
                       train.fraction = 0.5,       
                       mFeatures = 3,              
                       n.minobsinnode = 10,        
                       keep.data=TRUE,
                       cv.folds=1,                 
                       verbose = FALSE)            
  
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
  gbm.no.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,
                  data=data,               
                  var.monotone=c(0,0,0,0,0,0),
                  distribution= list(name="quantile",alpha=0.95),
                  n.trees=100,                 
                  shrinkage=0.005,             
                  interaction.depth=3,         
                  bag.fraction = 0.5,          
                  train.fraction = 0.5,        
                  mFeatures = 3,               
                  n.minobsinnode = 10,         
                  keep.data=TRUE,
                  cv.folds=1,                 
                  verbose = FALSE)            
  set.seed(15)
  gbm.zero.offset <- gbm(Y~X1+X2+X3+X4+X5+X6, 
                       data=data,             
                       var.monotone=c(0,0,0,0,0,0),
                       distribution= list(name="quantile",alpha=0.95), 
                       n.trees=100,                 
                       shrinkage=0.005,
                       interaction.depth=3,
                       bag.fraction = 0.5,
                       train.fraction = 0.5,
                       mFeatures = 3, 
                       n.minobsinnode = 10,
                       keep.data=TRUE,
                       offset=offset,
                       cv.folds=1,
                       verbose = FALSE)
  
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
  
  offset <- rep(0, N)
  
  # Generate new gbm object
  set.seed(15)
  gbm.no.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,    
                  data=data,                   
                  var.monotone=c(0,0,0,0,0,0), 
                  distribution="tdist",     
                  n.trees=100,              
                  shrinkage=0.005,          
                  interaction.depth=3,      
                  bag.fraction = 0.5,       
                  train.fraction = 0.5,     
                  mFeatures = 3,            
                  n.minobsinnode = 10,      
                  keep.data=TRUE,
                  cv.folds=1,               
                  verbose = FALSE)          
  set.seed(15)
  gbm.zero.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,
                       data=data,            
                       var.monotone=c(0,0,0,0,0,0),
                       distribution="tdist",     
                       n.trees=100,              
                       shrinkage=0.005,          
                       interaction.depth=3,      
                       bag.fraction = 0.5,       
                       train.fraction = 0.5,     
                       mFeatures = 3,            
                       n.minobsinnode = 10,      
                       keep.data=TRUE,
                       offset=offset,
                       cv.folds=1,               
                       verbose = FALSE)          
  
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
  gbm.no.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,
                  data=data,               
                  var.monotone=c(0,0,0,0,0,0),
                  distribution="tweedie",
                  n.trees=100,           
                  shrinkage=0.005,       
                  interaction.depth=3,   
                  bag.fraction = 0.5,    
                  train.fraction = 0.5,  
                  mFeatures = 3,         
                  n.minobsinnode = 10,   
                  keep.data=TRUE,
                  cv.folds=1,            
                  verbose = FALSE)       
  set.seed(15)
  gbm.zero.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,
                       data=data,            
                       var.monotone=c(0,0,0,0,0,0),
                       distribution="tweedie",    
                       n.trees=100,               
                       shrinkage=0.005,           
                       interaction.depth=3,       
                       bag.fraction = 0.5,        
                       train.fraction = 0.5,      
                       mFeatures = 3,             
                       n.minobsinnode = 10,       
                       keep.data=TRUE,
                       offset=offset,
                       cv.folds=1,                
                       verbose = FALSE)           
  
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
  gbm.smaller.offset <- gbm(Y~X1+X2+X3,        
                       data=data,              
                       weights=w,
                       offset=offset.smaller,
                       var.monotone=c(0,0,0),  
                       distribution="adaboost",
                       n.trees=100,            
                       shrinkage=0.001,        
                       interaction.depth=3,    
                       bag.fraction = 0.5,     
                       train.fraction = 0.5,   
                       cv.folds=1,             
                       n.minobsinnode = 10,    
                       verbose = FALSE)        
  set.seed(15)  
  gbm.larger.offset <- gbm(Y~X1+X2+X3,
                         data=data,                 
                         weights=w,
                         offset = offset.larger,
                         var.monotone=c(0,0,0),
                         distribution="adaboost",
                         n.trees=100, 
                         shrinkage=0.001,
                         interaction.depth=3,
                         bag.fraction = 0.5,        
                         train.fraction = 0.5,
                         cv.folds=1,                
                         n.minobsinnode = 10,       
                         verbose = FALSE)           
  
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
  gbm.smaller.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,
                            data=data,          
                            var.monotone=c(0,0,0,0,0,0),
                            distribution="poisson",
                            n.trees=100,        
                            offset=offset.smaller,
                            shrinkage=0.005,
                            interaction.depth=3,
                            bag.fraction = 0.5,
                            train.fraction = 0.5,
                            mFeatures = 3,
                            n.minobsinnode = 10,
                            keep.data=TRUE,
                            cv.folds=1,         
                            verbose = FALSE)    
  set.seed(15)
  gbm.larger.offset <- gbm(Y~X1+X2+X3+X4+X5+X6,
                         data=data,            
                         var.monotone=c(0,0,0,0,0,0),
                         distribution="Poisson",
                         n.trees=100,           
                         offset=offset.larger,
                         shrinkage=0.005,       
                         interaction.depth=3,   
                         bag.fraction = 0.5,    
                         train.fraction = 0.5,  
                         mFeatures = 3,         
                         n.minobsinnode = 10,   
                         keep.data=TRUE,
                         cv.folds=1,            
                         verbose = FALSE)       
    
  
  expect_true(gbm.smaller.offset$initF - gbm.larger.offset$initF > 0)
})
