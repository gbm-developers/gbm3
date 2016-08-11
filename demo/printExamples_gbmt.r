# Create some output to test the new print function and
# to be comparable with any future changes to gbmt.

library( MASS )

data( birthwt )
data( VA )
data( iris )
data( fgl )
data( cpus )
data( housing )

set.seed( 20090417 )

t_params <- training_params(num_trees = 1000, shrinkage = 0.01, num_features = ncol(birthwt)-1,
                            num_train = nrow(birthwt), id = seq_len(nrow(birthwt)))
bMod <- gbmt( low ~ ., data=birthwt, distribution=gbm_dist("Bernoulli"),
            train_params=t_params, cv_folds=5,
             is_verbose = FALSE)
bMod

bwt <- birthwt
bwt <- bwt[ sample( 1:nrow( bwt ) ),]
t_params <- training_params(num_trees = 1000, shrinkage = 0.01, num_features = ncol(bwt)-1,
                            num_train = round(0.9 * nrow(bwt)), id = seq_len(nrow(bwt)))
aMod <- gbmt( low ~ ., data=bwt, distribution=gbm_dist("AdaBoost"),
             train_params=t_params, cv_folds=10,
            is_verbose = FALSE )
aMod

t_params <- training_params(num_trees = 1000, shrinkage = 0.1, num_features = ncol(VA)-1,
                            num_train = nrow(VA), id = seq_len(nrow(VA)))
cMod <- gbmt( Surv( stime, status ) ~ treat + age + Karn + diag.time + cell + prior,
             data = VA, distribution=gbm_dist("CoxPH"), train_params=t_params, cv_folds = 5,
             is_verbose = FALSE)
cMod

t_params <- training_params(num_trees = 1000, shrinkage = 0.1, num_features = ncol(iris)-1,
                            num_train = round(0.9 * nrow(iris)), id = seq_len(nrow(iris)))
kMod <- gbmt( Species ~ . , data=iris , distribution = gbm_dist("Gaussian"), train_params=t_params,
             cv_folds=5)
kMod

t_params <- training_params(num_trees = 1000, shrinkage = 0.01, num_features = ncol(fgl)-1,
                            num_train = nrow(fgl), id = seq_len(nrow(fgl)))
kMod2 <- gbmt( type ~ ., data=fgl, distribution = gbm_dist("Gaussian"), train_params=t_params,
              cv_folds=5)
kMod2

mycpus <- cpus
mycpus <- mycpus[, -1 ]
t_params <- training_params(num_trees = 1000, shrinkage = 0.01, num_features = ncol(mycpus)-1,
                            num_train = nrow(mycpus), id = seq_len(nrow(mycpus)))
gMod <- gbmt( log( perf ) ~ ., data = mycpus, distribution=gbm_dist("Gaussian"),
             cv_folds=5, train_params = t_params,
             is_verbose = FALSE)
gMod

biMod <- gbmt( log(perf) ~ ., data=mycpus, distribution = gbm_dist("Gaussian"),
              cv_folds=5, train_params = t_params)
biMod

t_params <- training_params(num_trees = 1000, shrinkage = 0.01, num_features = ncol(mycpus)-1,
                            num_train = nrow(mycpus), id = seq_len(nrow(mycpus)),
                            interaction_depth = 3)
tMod <- gbmt( log(perf) ~ ., data=mycpus, distribution=gbm_dist("TDist"),
             cv_folds=5, train_params = t_params)
tMod

lMod <- gbmt( log(perf) ~ ., data=mycpus, distribution=gbm_dist("Laplace"),
             cv_folds=5, train_params = t_params)
lMod

qMod <- gbmt( log(perf) ~ ., data=mycpus,
             distribution=gbm_dist(name="Quantile", alpha=.7),
             cv_folds=5, train_params = t_params, is_verbose = FALSE)
qMod

t_params <- training_params(num_trees = 1000, shrinkage = 0.01, num_features = ncol(housing)-1,
                            num_train = nrow(housing), id = seq_len(nrow(housing)))
pMod <- gbmt( Freq ~ ., data=housing , distribution=gbm_dist("Poisson"),
              cv_folds=5 , train_params = t_params)
pMod