
# COX PROPORTIONAL HAZARDS REGRESSION EXAMPLE

cat("Running cox proportional hazards regression example.\n")
# create some data
N <- 3000
X1 <- runif(N)
X2 <- runif(N)
X3 <- factor(sample(letters[1:4],N,replace=T))
mu <- c(-1,0,1,2)[as.numeric(X3)]

f <- 0.5*sin(3*X1 + 5*X2^2 + mu/10)
tt.surv <- rexp(N,exp(f))
tt.cens <- rexp(N,0.5)
delta <- as.numeric(tt.surv <= tt.cens)
tt <- apply(cbind(tt.surv,tt.cens),1,min)

# throw in some missing values
X1[sample(1:N,size=100)] <- NA
X3[sample(1:N,size=300)] <- NA

# random weights if you want to experiment with them
w <- rep(1,N)

data <- data.frame(tt=tt,delta=delta,X1=X1,X2=X2,X3=X3)

# fit initial model
gbm1 <- gbm(Surv(tt,delta)~X1+X2+X3,       # formula
            data=data,                 # dataset
            weights=w,
            var.monotone=c(0,0,0),     # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
            distribution="coxph",
            n.trees=3000,              # number of trees
            shrinkage=0.001,           # shrinkage or learning rate, 0.001 to 0.1 usually work
            interaction.depth=3,       # 1: additive model, 2: two-way interactions, etc
            bag.fraction = 0.5,        # subsampling fraction, 0.5 is probably best
            train.fraction = 0.5,      # fraction of data for training, first train.fraction*N used for training
            cv.folds = 5,              # do 5-fold cross-validation
            n.minobsinnode = 10,       # minimum total weight needed in each node
            keep.data = TRUE,
            verbose = FALSE)           # don't print progress

# plot the performance
best.iter <- gbm.perf(gbm1,method="OOB")  # returns out-of-bag estimated best number of trees
print(best.iter)
best.iter <- gbm.perf(gbm1,method="cv") # returns test set estimate of best number of trees
print(best.iter)
best.iter <- gbm.perf(gbm1,method="test") # returns test set estimate of best number of trees
print(best.iter)

# plot variable influence
summary(gbm1,n.trees=1)         # based on the first tree
summary(gbm1,n.trees=best.iter) # based on the estimated best number of trees

# create marginal plots
# plot variable X1,X2,X3 after "best" iterations
par(mfrow=c(1,3))
plot.gbm(gbm1,1,best.iter)
plot.gbm(gbm1,2,best.iter)
plot.gbm(gbm1,3,best.iter)
par(mfrow=c(1,1))
plot.gbm(gbm1,1:2,best.iter) # contour plot of variables 1 and 2 after "best" number iterations

# 3-way plots
plot.gbm(gbm1,1:3,best.iter)

# print the first and last trees... just for curiosity
pretty.gbm.tree(gbm1,1)
pretty.gbm.tree(gbm1,gbm1$n.trees)

# make some new data
N <- 1000
X1 <- runif(N)
X2 <- runif(N)
X3 <- factor(sample(letters[1:4],N,replace=T))
mu <- c(-1,0,1,2)[as.numeric(X3)]

f <- 0.5*sin(3*X1 + 5*X2^2 + mu/10)
tt.surv <- rexp(N,exp(f))
tt.cens <- rexp(N,0.5)

data2 <- data.frame(tt=apply(cbind(tt.surv,tt.cens),1,min),
                    delta=as.numeric(tt.surv <= tt.cens),
                    X1=X1,X2=X2,X3=X3)

# predict on the new data using "best" number of trees
# f.predict will be on the canonical scale (logit,log,etc.)
f.predict <- predict(gbm1,data2,best.iter)

# Cox PH error
# boosting
risk <- rep(0,N)
for(i in 1:N)
{
   risk[i] <- sum( (data2$tt>=data2$tt[i])*exp(f.predict) )
}
cat("Boosting:",sum( data2$delta*( f.predict - log(risk) ) ),"\n")

# linear model
coxph1 <- coxph(Surv(tt,delta)~X1+X2+X3,data=data)
f.predict <- predict(coxph1,newdata=data2)
risk <- rep(0,N)
for(i in 1:N)
{
   risk[i] <- sum( (data2$tt>=data2$tt[i])*exp(f.predict) )
}
cat("Linear model:",sum( data2$delta*( f.predict - log(risk) ) ),"\n")
