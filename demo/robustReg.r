library(MASS)

set.seed(20090415)


x <- mvrnorm(100, mu=rep(0,5), Sigma=diag(rep(1,5)))
r <- rnorm(100)
r <- ifelse(runif(100) < 0.25, r*4, r)
y <- apply(x, 1, sum) + r

d <- data.frame(y=y, x)

gmod <- gbm(y ~ ., data=d, distribution="gaussian",
            n.tree = 2000, shrinkage = 0.01, cv.folds=5,
            verbose = FALSE)
tmod4 <- gbm(y ~ ., data=d, distribution="tdist", # defaults to 4 df
             n.tree=2000, shrinkage = 0.01, cv.folds=5,
             verbose = FALSE)
tmod6 <- gbm(y ~ ., data=d, distribution=list(name="tdist", df=6),
             n.tree=2000, shrinkage = 0.01, cv.folds=5,
             verbose = FALSE)
tmod100 <- gbm(y ~ ., data=d, distribution=list(name="tdist", df=100),
               n.tree=2000, shrinkage = 0.01, cv.folds=5,
               verbose = FALSE)

par.old <- par(mfrow=c(2,2))
gbest <- gbmt_performance(gmod, method="cv")
t4best <- gbmt_performance(tmod4, method="cv")
t6best <- gbmt_performance(tmod6, method="cv")
t100best <- gbmt_performance(tmod100, method="cv")
par(par.old)