library( MASS )

set.seed( 20090415 )


x <- mvrnorm( 100, mu=rep( 0, 5 ) , Sigma=diag( rep( 1, 5 ) ) )
r <- rnorm( 100 )
r <- ifelse( runif( 100 ) < .25 , r * 4, r )
y <- apply( x, 1, sum ) + r

d <- data.frame( y=y , x)
train_params <- training_params(num_trees = 2000, interaction_depth = 1, bag_fraction = 0.5,
                                num_train = nrow(d), id = seq_len(nrow(d)), num_features = ncol(x),
                                shrinkage = 0.01)

gmod <- gbmt( y ~ ., data=d, distribution=gbm_dist("Gaussian"),
             train_params = train_params, cv_folds=5, is_verbose = FALSE)
tmod4 <- gbmt( y ~ ., data=d, distribution=gbm_dist("TDist"), # defaults to 4 df
             train_params = train_params, cv_folds=5, is_verbose = FALSE)
tmod6 <- gbmt( y ~ ., data=d, distribution=gbm_dist( name="TDist", df=6 ),
              train_params = train_params,cv_folds=5, is_verbose = FALSE)
tmod100 <- gbmt( y ~ ., data=d, distribution=gbm_dist( name="TDist", df=100 ),
               train_params = train_params, cv_folds=5, is_verbose = FALSE)

par(mfrow=c( 2, 2 ) )
gbest <- gbm_perf( gmod , method="cv" )
t4best <- gbm_perf( tmod4 , method="cv" )
t6best <- gbm_perf( tmod6 , method="cv" )
t100best <- gbm_perf( tmod100 , method="cv" )

dev.off()