context("testing openmp parallelization")
test_that("gbm refuses to work with insane numbers of threads", {
  skip_on_os("mac")
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
  Y <- Y + rnorm(N, 0, sigma)
  
  ## create a bunch of missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  w <- rep(1,N)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # fit initial model
  expect_error(gbmt(Y~X1+X2+X3+X4+X5+X6,         # formula
                   data=data,                   # dataset
                   var_monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                   keep_gbm_data=TRUE,
                   cv_folds=10, # do 10-fold cross-validation
                   par_details=gbmParallel(num_threads=-1)),
               "number of threads must be strictly positive",
               fixed=TRUE)
})
test_that("gbm refuses to work with insane array chunk size - old api", {
  skip_on_os("mac")
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
  
  ## create a bunch of missing values
  X1[sample(1:N,size=100)] <- NA
  X3[sample(1:N,size=300)] <- NA
  
  w <- rep(1,N)
  
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # fit initial model
  expect_error(gbmt(Y~X1+X2+X3+X4+X5+X6,         # formula
                   data=data,                   # dataset
                   var_monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                   keep_gbm_data=TRUE,
                   cv_folds=10, # do 10-fold cross-validation
                   par_details=gbmParallel(num_threads=1, array_chunk_size=0)),
               "array chunk size must be strictly positive", fixed=TRUE)
})
test_that("gbm refuses to work with insane numbers of threads - old API", {
  skip_on_os("mac")
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

    ## create a bunch of missing values
    X1[sample(1:N,size=100)] <- NA
    X3[sample(1:N,size=300)] <- NA

    w <- rep(1,N)

    data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)

    # fit initial model
    expect_error(gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
                     data=data,                   # dataset
                     var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                     distribution="Gaussian",     # bernoulli, adaboost, gaussian, poisson, coxph, or
                                        # list(name="quantile",alpha=0.05) for quantile regression
                     weights=rep(1,nrow(data)),
                     n.trees=2000,                 # number of trees
                     shrinkage=0.005,             # shrinkage or learning rate, 0.001 to 0.1 usually work
                     interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc.
                     bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
                     train.fraction = 0.5,        # fraction of data for training, first train.fraction*N used for training
                     n.minobsinnode = 10,         # minimum number of obs needed in each node
                     keep.data=TRUE,
                     cv.folds=10, # do 10-fold cross-validation
                     par.details=gbmParallel(num_threads=-1)),
                 "number of threads must be strictly positive",
                 fixed=TRUE)
})
test_that("gbm refuses to work with insane array chunk size - old api", {
  skip_on_os("mac")
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

    ## create a bunch of missing values
    X1[sample(1:N,size=100)] <- NA
    X3[sample(1:N,size=300)] <- NA

    w <- rep(1,N)

    data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)

    # fit initial model
    expect_error(gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
                     data=data,                   # dataset
                     var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                     distribution="Gaussian",     # bernoulli, adaboost, gaussian, poisson, coxph, or
                                        # list(name="quantile",alpha=0.05) for quantile regression
                     n.trees=2000,                 # number of trees
                     shrinkage=0.005,             # shrinkage or learning rate, 0.001 to 0.1 usually work
                     interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc.
                     bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
                     train.fraction = 0.5,        # fraction of data for training, first train.fraction*N used for training
                     n.minobsinnode = 10,         # minimum number of obs needed in each node
                     keep.data=TRUE,
                     cv.folds=10, # do 10-fold cross-validation
                     par.details=gbmParallel(num_threads=2, array_chunk_size=0)),
                 "array chunk size must be strictly positive", fixed=TRUE)
})
