library( MASS )

set.seed( 20090415 )


x <- mvrnorm( 100, mu=rep( 0, 5 ) , Sigma=diag( rep( 1, 5 ) ) )
r <- rnorm( 100 )
r <- ifelse( runif( 100 ) < .25 , r * 4, r )
y <- apply( x, 1, sum ) + r

d <- data.frame( y=y , x)

gmod <- gbm( y ~ ., data=d, distribution="Gaussian",
             n.tree = 2000, shrinkage = .01 , cv.folds=5,
            verbose = FALSE)
tmod4 <- gbm( y ~ ., data=d, distribution="TDist", # defaults to 4 df
              n.tree=2000, shrinkage = .01, cv.folds=5,
             verbose = FALSE)
tmod6 <- gbm( y ~ ., data=d, distribution=list( name="TDist", df=6 ),
              n.tree=2000, shrinkage = .01, cv.folds=5,
              verbose = FALSE)
tmod100 <- gbm( y ~ ., data=d, distribution=list( name="TDist", df=100 ),
              n.tree=2000, shrinkage = .01, cv.folds=5,
               verbose = FALSE)

par(mfrow=c( 2, 2 ) )
gbest <- gbm_perf( gmod , method="cv" )
t4best <- gbm_perf( tmod4 , method="cv" )
t6best <- gbm_perf( tmod6 , method="cv" )
t100best <- gbm_perf( tmod100 , method="cv" )

dev.off()
