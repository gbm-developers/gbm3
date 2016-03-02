context("fold.id")

test_that("gaussian works in parallel", {
    skip_on_cran()
    ## Based on example in R package
    # Data generation from gbm manual ------------------------------------------
    set.seed(1)
    N     <- 1000
    X1    <- runif(N)
    X2    <- 2*runif(N)
    X3    <- ordered(sample(letters[1:4],N,replace=TRUE),levels=letters[4:1])
    X4    <- factor(sample(letters[1:6],N,replace=TRUE))
    X5    <- factor(sample(letters[1:3],N,replace=TRUE))
    X6    <- 3*runif(N)
    mu    <- c(-1,0,1,2)[as.numeric(X3)]
    SNR   <- 10
    Y     <- X1**1.5 + 2 * (X2**.5) + mu
    sigma <- sqrt(var(Y)/SNR)
    Y     <- Y + rnorm(N,0,sigma)
    X1[sample(1:N,size=500)] <- NA
    X4[sample(1:N,size=300)] <- NA
    data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
     
    # Fit some models with fold.id ---------------------------------------------
    set.seed(1)
    folds <- sample(seq(1,3), nrow(data), replace = TRUE)
     
    set.seed(123456)
    gid1 <- gbm(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                data              = data,
                distribution      = "gaussian",
                n.trees           = 5,
                bag.fraction      = 1.0,
                shrinkage         = 0.01,
                interaction.depth = 3,
                fold.id           = folds)
     
    set.seed(260591)
    gid2 <- gbm(Y ~ X1 + X2 + X3 + X4 + X5 + X6,
                data              = data,
                distribution      = "gaussian",
                n.trees           = 5,
                bag.fraction      = 1.0,
                shrinkage         = 0.01,
                interaction.depth = 3,
                fold.id           = folds)
     

    expect_true(all(gid1$fold.id == gid2$fold.id))
})