# LOGISTIC REGRESSION EXAMPLE

cat("Running logistic regression example.\n")
# create some data
N <- 1000
X1 <- runif(N)
X2 <- runif(N)
X3 <- factor(sample(letters[1:4],N,replace=TRUE))
mu <- c(-1,0,1,2)[as.numeric(X3)]

p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
Y <- rbinom(N,1,p)

# random weights if you want to experiment with them
w <- rexp(N)
w <- N*w/sum(w)

data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)

train_params <- training_params(num_trees = 3000, 
                                shrinkage = 0.001, 
                                bag_fraction = 0.5,
                                num_train = N/2, 
                                id=seq_len(nrow(data)), 
                                min_num_obs_in_node = 10,
                                interaction_depth = 3, 
                                num_features = 3)

# fit initial model
gbm1 <- gbmt(Y~X1+X2+X3,                # formula
            data=data,                  # dataset
            weights=w,
            var_monotone=c(0,0,0),      # -1: monotone decrease, 
                                        # +1: monotone increase, 
                                        # 0: no monotone restrictions
            distribution=gbm_dist("Bernoulli"),
            train_params = train_params,
            cv_folds=5,                 # do 5-fold cross-validation
            is_verbose = FALSE)         # don't print progress

# plot the performance
# returns out-of-bag estimated best number of trees
best.iter.oob <- gbmt_performance(gbm1,method="OOB")  
print(best.iter.oob)
# returns 5-fold cv estimate of best number of trees
best.iter.cv <- gbmt_performance(gbm1,method="cv")   
print(best.iter.cv)
# returns test set estimate of best number of trees
best.iter.test <- gbmt_performance(gbm1,method="test") 
print(best.iter.test)

best.iter <- best.iter.test

# plot variable influence
summary(gbm1,num_trees=1)         # based on the first tree
summary(gbm1,num_trees=best.iter) # based on the estimated best number of trees

# create marginal plots
# plot variable X1,X2,X3 after "best" iterations
par.old <- par(mfrow=c(1,3))
plot(gbm1,1,best.iter)
plot(gbm1,2,best.iter)
plot(gbm1,3,best.iter)
par(mfrow=c(1,1))
plot(gbm1,1:2,best.iter) # contour plot of variables 1 and 2 after "best" number iterations
plot(gbm1,2:3,best.iter) # lattice plot of variables 2 and 3 after "best" number iterations

# 3-way plot
plot(gbm1,1:3,best.iter)

# print the first and last trees
print(pretty_gbm_tree(gbm1,1))
print(pretty_gbm_tree(gbm1, gbm1$params$num_trees))

# make some new data
N <- 1000
X1 <- runif(N)
X2 <- runif(N)
X3 <- factor(sample(letters[1:4],N,replace=TRUE))
mu <- c(-1,0,1,2)[as.numeric(X3)]

p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
Y <- rbinom(N,1,p)
data2 <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)

# predict on the new data using "best" number of trees
# f.predict will be on the canonical scale (logit,log,etc.)
f.predict <- predict(gbm1,data2,
                     num_trees=c(best.iter.oob,best.iter.cv,best.iter.test))
# transform to probability scale for logistic regression
p.pred <- 1/(1+exp(-f.predict))

# calibration plot for logistic regression - well calibrated means a 45 degree line
par(mfrow=c(1,1))
calibrate_plot(Y,p.pred[,3])
par(par.old) # reset graphics options to previous settings

# logistic error
sum(data2$Y*f.predict[,1] - log(1+exp(f.predict[,1])))
sum(data2$Y*f.predict[,2] - log(1+exp(f.predict[,2])))
sum(data2$Y*f.predict[,3] - log(1+exp(f.predict[,3])))
