# RANKING EXAMPLE

cat("Running ranking (LambdaMart) example.\n")

# Create synthetic data that shows how pairwise training can be better
# Note: no claim to represent 'real world' data!

generate.data <- function(N) {
  
  # create query groups, with an average size of 25 items each
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
  
  data.frame(Y, query=query, X1, X2, X3)
}

cat('Generating data\n')
N=1000
data.train <- generate.data(N)

# Now we fit 3 different models to the same data:
# * Gaussian
# * Pairwise with NDCG ranking metric
# * Pairwise with CONC (fraction of concordant pairs) ranking metric

cat('Fitting a model with gaussian loss function\n')
train_params_gauss <- training_params(num_trees = 2000, shrinkage = 0.005, 
                                      interaction_depth = 3, bag_fraction = 0.5,
                                      num_train = N, id = seq_len(nrow(data)), num_features = 3,
                                      min_num_obs_in_node = 10)
gbm.gaussian <- gbmt(Y~X1+X2+X3,      # formula
                    data=data.train,         # dataset
                    distribution=gbm_dist('Gaussian'), # loss function: gaussian
                    train_params=train_params_gauss,
                    keep_gbm_data=TRUE,          # store copy of input data in model
                    cv_folds=5,              # number of cross validation folds
                    is_verbose = FALSE         # don't print progress
)             # use a single core (to prevent possible problems caused by wronly detecting cores)

# estimate number of trees
best.iter.gaussian <- gbm_perf(gbm.gaussian, method="cv")
title('Training of gaussian model')

cat('Fitting a model with pairwise loss function (ranking metric: normalized discounted cumulative gain)\n')

gbm.ndcg <- gbmt(Y~X1+X2+X3,          # formula
                data=data.train,     # dataset
                distribution=gbm_dist(  # loss function:
                  name='Pairwise',   # pairwise
                  metric="ndcg",     # ranking metric: normalized discounted cumulative gain
                  group='query'),    # column indicating query groups
                train_params=train_params_gauss,
                keep_gbm_data=TRUE,      # store copy of input data in model
                cv_folds=5,          # number of cross validation folds
                is_verbose = FALSE     # don't print progress
)         # use a single core
# estimate number of trees
best.iter.ndcg <- gbm_perf(gbm.ndcg, method='cv')
title('Training of pairwise model with ndcg metric')

cat('Fit a model with pairwise loss function (ranking metric: fraction of concordant pairs)\n')

gbm.conc <- gbmt(Y~X1+X2+X3,          # formula
                data=data.train,     # dataset
                distribution=gbm_dist(   # loss function:
                  name='Pairwise',   # pairwise
                  metric="conc",     # ranking metric: concordant pairs
                  group='query'),    # column indicating query groups
                train_params = train_params_gauss,
                keep_gbm_data=TRUE,      # store copy of input data in model
                cv_folds=5,          # number of cross validation folds
                is_verbose = FALSE     # don't print progress
)         # use a single core

# estimate number of trees
best.iter.conc <- gbm_perf(gbm.conc, method='cv')
title('Training of pairwise model with conc metric')


# plot variable importance
par.old <- par(mfrow=c(1,3))
summary(gbm.gaussian, num_trees=best.iter.gaussian, main='gaussian')
summary(gbm.ndcg, num_trees=best.iter.ndcg, main='pairwise (ndcg)')
summary(gbm.conc, num_trees=best.iter.conc, main='pairwise (conc)')
par(par.old)

cat("Generating some new data\n")

data.test <- generate.data(N)

cat("Calculating predictions\n")

predictions <- data.frame(random=runif(N),
                          X2=data.test$X2,
                          gaussian=predict(gbm.gaussian, data.test, best.iter.gaussian),
                          pairwise.ndcg=predict(gbm.ndcg, data.test, best.iter.ndcg),
                          pairwise.conc=predict(gbm.conc, data.test, best.iter.conc))

cat("Computing loss metrics\n")
dist_1 <- gbm_dist("Gaussian")
dist_2 <- gbm_dist(name='Pairwise', metric="ndcg", 
                   group_index=data.test$query, max_rank=5)
dist_3 <- gbm_dist(name='Pairwise', metric="conc",
                   group_index=data.test$query, max_rank=0)

result.table <- data.frame(measure=c('random', 'X2 only', 'gaussian', 'pairwise (ndcg)', 'pairwise (conc)'),
                           squared.loss=sapply(1:length(predictions), FUN=function(i) {
                             loss(y=data.test$Y, predictions[[i]], w=rep(1,N), offset=rep(0, N), dist=gbm_dist(name="Gaussian")) }),
                           ndcg5.loss=sapply(1:length(predictions), FUN=function(i) {
                             loss(y=data.test$Y, predictions[[i]], w=rep(1,N), offset=rep(0 ,N), distribution_obj=dist_2)}),
                           concordant.pairs.loss=sapply(1:length(predictions), FUN=function(i) {
                             loss(y=data.test$Y, predictions[[i]], w=rep(1,N), offset=rep(0, N), distribution_obj= dist_3) }),
                           row.names=NULL)

cat('Performance measures for the different models on the test set (smaller is better):\n')
print(result.table,digits=2)

# Brief explanation: Variable X1 is not correlated with the order of items, only
# with queries. Variable X2 is the only one that is correlated with the order of
# items within queries. However, it has a high query-correlated variance.
# Therefore, the 'optimal' possible ranking is just by X2. Of course, the
# pairwise models don't know this and don't completely achieve the same
# accuracy, due to noise and data limitation.
#
# The Gaussian model uses mostly X1, due to the high variance of X2; on the
# contrary, the pairwise models rely mainly on X2. The loss table shows that
# both pairwise models are better in terms of the ranking metrics, but worse in
# terms of squared loss.
