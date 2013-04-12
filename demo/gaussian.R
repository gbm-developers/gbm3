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

# fit initial model
gbm1 <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
            data=data,                   # dataset
            var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
            distribution="gaussian",     # bernoulli, adaboost, gaussian, poisson, coxph, or
                                         # list(name="quantile",alpha=0.05) for quantile regression
            n.trees=2000,                 # number of trees
            shrinkage=0.005,             # shrinkage or learning rate, 0.001 to 0.1 usually work
            interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc
            bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
            train.fraction = 0.5,        # fraction of data for training, first train.fraction*N used for training
            n.minobsinnode = 10,         # minimum number of obs needed in each node
            keep.data=TRUE,
            cv.folds=10,                 # do 10-fold cross-validation
            verbose = FALSE)             # don't print progress

# plot the performance
best.iter <- gbm.perf(gbm1,method="OOB")  # returns out-of-bag estimated best number of trees
best.iter <- gbm.perf(gbm1,method="test") # returns test set estimate of best number of trees
best.iter <- gbm.perf(gbm1,method="cv")   # returns cv estimate of best number of trees

# plot variable influence
summary(gbm1,n.trees=1)         # based on the first tree
summary(gbm1,n.trees=best.iter) # based on the estimated best number of trees

# print the first and last trees
print(pretty.gbm.tree(gbm1,1))
print(pretty.gbm.tree(gbm1,gbm1$n.trees))

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
plot(gbm1,1,best.iter)
plot(gbm1,2,best.iter)
plot(gbm1,3,best.iter)
par(mfrow=c(1,1))
plot(gbm1,1:2,best.iter) # contour plot of variables 1 and 2 after "best" number iterations
plot(gbm1,2:3,best.iter) # lattice plot of variables 2 and 3 after "best" number iterations
plot(gbm1,3:4,best.iter) # lattice plot of variables 2 and 3 after "best" number iterations

plot(gbm1,c(1,2,6),best.iter,cont=20) # 3-way plots
plot(gbm1,1:3,best.iter)
plot(gbm1,2:4,best.iter)
plot(gbm1,3:5,best.iter)

# check interactions
interact.gbm(gbm1,data=data,i.var=1:2,n.trees=best.iter)
# get all two way interactions
i.var <- subset(expand.grid(x1=1:6,x2=1:6), x1<x2)
rownames(i.var) <- apply(i.var,1,paste,collapse=":",sep="")
apply(i.var,1,
      function(i.var) interact.gbm(gbm1,data=data,i.var=i.var,n.trees=best.iter))
