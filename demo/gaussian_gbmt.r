# LEAST SQUARES EXAMPLE

cat("Running least squares regression example.\n")
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

# random weights if you want to experiment with them
# w <- rexp(N)
# w <- N*w/sum(w)
w <- rep(1,N)

data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
require(gbm)
train_params <- training_params(num_trees = 2000, bag_fraction = 0.5, num_train = N/2,
                                id = seq_len(nrow(data)),
                                min_num_obs_in_node = 10, num_features = 3, 
                                interaction_depth = 3)
# fit initial model
gbm1 <- gbmt(Y~X1+X2+X3+X4+X5+X6,         # formula
            data=data,                   # dataset
            var_monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
            distribution=gbm_dist("Gaussian"),     # bernoulli, adaboost, gaussian, poisson, coxph, or
            # list(name="quantile",alpha=0.05) for quantile regression
            train_params=train_params,
            keep_gbm_data=TRUE,
            cv_folds=10,                 # do 10-fold cross-validation
            is_verbose = FALSE)             # don't print progress

str(gbm1,max.level=1)
# plot the performance
best.iter <- gbm_perf(gbm1,method="OOB")  # returns out-of-bag estimated best number of trees
best.iter <- gbm_perf(gbm1,method="test") # returns test set estimate of best number of trees
best.iter <- gbm_perf(gbm1,method="cv")   # returns cv estimate of best number of trees

# plot variable influence
summary(gbm1,num_trees=1)         # based on the first tree
summary(gbm1,num_trees=best.iter) # based on the estimated best number of trees

# print the first and last trees
print(pretty_gbm_tree(gbm1,1))
print(pretty_gbm_tree(gbm1,gbm1$params$num_trees))

print(gbm1$c.splits[1:3])

# make some new data
N <- 1000
X1 <- runif(N)
X2 <- 2*runif(N)
X3 <- factor(sample(letters[1:4],N,replace=TRUE))
X4 <- ordered(sample(letters[1:6],N,replace=TRUE))
X5 <- factor(sample(letters[1:3],N,replace=TRUE))
X6 <- 3*runif(N)
mu <- c(-1,0,1,2)[as.numeric(X3)]

Y <- X1**1.5 + 2 * (X2**.5) + mu
Y <- Y + rnorm(N,0,sigma)

data2 <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
print(data2[1:10,])

# predict on the new data using "best" number of trees
f.predict <- predict(gbm1,data2,best.iter) # f.predict will be on the canonical scale (logit,log,etc.)

print(f.predict[1:10])
# least squares error
print(sum((data2$Y-f.predict)^2))

# create marginal plots
# plot variable X1,X2,X3 after "best" iterations
par(mfrow=c(1,3))
gbmt_plot(gbm1,1,best.iter)
gbmt_plot(gbm1,2,best.iter)
gbmt_plot(gbm1,3,best.iter)
par(mfrow=c(1,1))
gbmt_plot(gbm1,1:2,best.iter) # contour plot of variables 1 and 2 after "best" number iterations
gbmt_plot(gbm1,2:3,best.iter) # lattice plot of variables 2 and 3 after "best" number iterations
gbmt_plot(gbm1,3:4,best.iter) # lattice plot of variables 2 and 3 after "best" number iterations

gbmt_plot(gbm1,c(1,2,6),best.iter,cont=20) # 3-way plots
gbmt_plot(gbm1,1:3,best.iter)
gbmt_plot(gbm1,2:4,best.iter)
gbmt_plot(gbm1,3:5,best.iter)

# check interactions
interact(gbm1,data=data,var_indices=1:2, num_trees=best.iter)
# get all two way interactions
i.var <- subset(expand.grid(x1=1:6,x2=1:6), x1<x2)
rownames(i.var) <- apply(i.var,1,paste,collapse=":",sep="")
apply(i.var,1,
      function(i.var) interact(gbm1,data=data,var_indices=i.var,num_trees=best.iter))