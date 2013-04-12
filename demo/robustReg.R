library( MASS )

set.seed( 20090415 )


x <- mvrnorm( 100, mu=rep( 0, 5 ) , Sigma=diag( rep( 1, 5 ) ) )
r <- rnorm( 100 )
r <- ifelse( runif( 100 ) < .25 , r * 4, r )
y <- apply( x, 1, sum ) + r

d <- data.frame( y=y , x)

gmod <- gbm( y ~ ., data=d, distribution="gaussian",
             n.tree = 2000, shrinkage = .01 , cv.folds=5,
            verbose = FALSE, n.cores=1)
tmod4 <- gbm( y ~ ., data=d, distribution="tdist", # defaults to 4 df
              n.tree=2000, shrinkage = .01, cv.folds=5,
             verbose = FALSE, n.cores=1)
tmod6 <- gbm( y ~ ., data=d, distribution=list( name="tdist", df=6 ),
              n.tree=2000, shrinkage = .01, cv.folds=5,
              verbose = FALSE, n.cores=1)
tmod100 <- gbm( y ~ ., data=d, distribution=list( name="tdist", df=100 ),
              n.tree=2000, shrinkage = .01, cv.folds=5,
               verbose = FALSE, n.cores=1)

par(mfrow=c( 2, 2 ) )
gbest <- gbm.perf( gmod , method="cv" )
t4best <- gbm.perf( tmod4 , method="cv" )
t6best <- gbm.perf( tmod6 , method="cv" )
t100best <- gbm.perf( tmod100 , method="cv" )

qscale <- function( x ){
  x / abs( diff( quantile( x , prob=c( .25, .75 ) ) ) )
}

rg <- qscale( resid( gmod , n.trees=gbest) )
rt4 <- qscale( resid( tmod4 , n.trees=t4best) )
rt6 <- qscale( resid( tmod6 , n.trees=t6best) )
rt100 <- qscale( resid( tmod100 , n.trees=t100best ) )

ylimits <- range(rg, rt4, rt6, rt100)

plot( rg, main="Gaussian", ylim=ylimits ); abline( h=0 )
plot( rt4, main="t(4)", ylim=ylimits ); abline( h=0 )
plot( rt6, main="t(6)", ylim=ylimits ); abline( h=0 )
plot( rt100, main="t(100)", ylim=ylimits ); abline( h=0 )

dev.off()
