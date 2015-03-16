context("test parallelization")

test_that("gaussian works in parallel", {
    skip_on_cran()
    ## Based on example in R package
    
    ## test Gaussian distribution gbm model
    set.seed(1)

    # create some data
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

    # fit initial model
    gbm1 <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
                data=data,                   # dataset
                var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                distribution="gaussian",     # bernoulli, adaboost, gaussian, poisson, coxph, or
                                            # list(name="quantile",alpha=0.05) for quantile regression
                n.trees=2000,                 # number of trees
                shrinkage=0.005,             # shrinkage or learning rate, 0.001 to 0.1 usually work
                interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc.
                bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
                train.fraction = 0.5,        # fraction of data for training, first train.fraction*N used for training
                n.minobsinnode = 10,         # minimum number of obs needed in each node
                keep.data=TRUE,
                cv.folds=10, # do 10-fold cross-validation
                n.cores=2)                 

    # Get best model
    best.iter <- gbm.perf(gbm1,method="cv")   # returns cv estimate of best number of trees
    
    set.seed(2)
    # make some new data
    N <- 1000
    X1 <- runif(N)
    X2 <- 2*runif(N)
    X3 <- factor(sample(letters[1:4],N,replace=TRUE))
    X4 <- ordered(sample(letters[1:6],N,replace=TRUE))
    X5 <- factor(sample(letters[1:3],N,replace=TRUE))
    X6 <- 3*runif(N)
    mu <- c(-1,0,1,2)[as.numeric(X3)]

    # Actual underlying signal
    Y <- X1**1.5 + 2 * (X2**.5) + mu

    # Want to see how close predictions are to the underlying signal; noise would just interfere with this
    # Y <- Y + rnorm(N,0,sigma)
    data2 <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)

    # predict on the new data using "best" number of trees
    f.predict <- predict(gbm1,data2,best.iter) # f.predict will be on the canonical scale (logit,log,etc.)

    # Base the validation tests on observed discrepancies
    expect_true(cor(data2$Y, f.predict) > 0.990)
    expect_true(sd(data2$Y-f.predict) < sigma)
})
