####################
# Author: James Hickey
#
# Series of tests to check that the basic gbm behaviour
#
####################

context("some basic checks")
test_that("gaussian works", {
  ## Based on example in R package
  
  ## test Gaussian distribution gbm model
  set.seed(1)
  
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
  
  w <- rep(1,N)
  offset <- rep(0, N)
  data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  
  # Set up for new API
  params <- training_params(num_trees=2000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.005, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=6)
  dist <- gbm_dist("Gaussian")
  
  fit <- gbmt(Y~X1+X2+X3+X4+X5+X6, data=data, distribution=dist, weights=w, offset=offset,
              train_params=params, var_monotone=c(0, 0, 0, 0, 0, 0), keep_gbm_data=TRUE, cv_folds=10, is_verbose=FALSE)
  
  best_iter <- gbm_perf(fit, method="cv") # returns test set estimate of best number of trees
  
  # Make prediction
  set.seed(2)
  # make some new data
  N <- 1000
  X1 <- runif(N)
  X2 <- 2*runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=TRUE))
  X4 <- ordered(sample(letters[1:6],N,replace=TRUE))
  X5 <- factor(sample(letters[1:3],N,replace=TRUE))
  X6 <- 3*runif(N)
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  # Actual underlying signal - Check how close to this we are
  Y <- X1**1.5 + 2 * (X2**.5) + mu
  data2 <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
  
  # predict on the new data using "best" number of trees
  f.predict <- predict(fit, data2, best_iter) # f.predict will be on the canonical scale (logit,log,etc.)
  
  # Base the validation tests on observed discrepancies
  expect_true(cor(data2$Y, f.predict) > 0.990)
  expect_true(sd(data2$Y-f.predict) < sigma)
})
test_that("coxph works - efron", {
  # Require Surv to be available
  require(survival)
  
  # create some data
  set.seed(1)
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
  
  # Put into new API
  dist <- gbm_dist("CoxPH", prior_node_coeff_var=10)
  params <- training_params(num_trees=3000, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=3)
  
  # fit initial model
  gbm1 <- gbmt(Surv(tt,delta)~X1+X2+X3, data=data, distribution=dist, weights=w, offset=rep(0, N),
               train_params=params, var_monotone=c(0, 0, 0), keep_gbm_data=TRUE, cv_folds=5, is_verbose=FALSE)
  
  best_iter <- gbm_perf(gbm1, method="test") # returns test set estimate of best number of trees
  
  # make some new data
  set.seed(2)
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4],N,replace=T))
  mu <- c(-1,0,1,2)[as.numeric(X3)]
  
  f <- 0.5*sin(3*X1 + 5*X2^2 + mu/10)  # -0.5 <= f <= 0.5 via sin fn.
  tt.surv <- rexp(N,exp(f))
  tt.cens <- rexp(N,0.5)
  
  data2 <- data.frame(tt=apply(cbind(tt.surv,tt.cens),1,min),
                      delta=as.numeric(tt.surv <= tt.cens),
                      f=f,
                      X1=X1,X2=X2,X3=X3)
  
  # predict on the new data using "best" number of trees
  # f.predict will be on the canonical scale (logit,log,etc.)
  f.predict <- predict(gbm1, data2, best_iter)
  
  #plot(data2$f,f.predict)
  # Use observed sd
  expect_true(sd(data2$f - f.predict) < 0.4)
})
test_that("coxph works - breslow", {
    # Require Surv to be available
    require(survival)
  
    # create some data
    set.seed(1)
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

    # Put into new API
    dist <- gbm_dist("CoxPH", prior_node_coeff_var=10, ties="breslow")
    params <- training_params(num_trees=3000, interaction_depth=3, min_num_obs_in_node=10, 
                              shrinkage=0.001, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=3)
    
    # fit initial model
    gbm1 <- gbmt(Surv(tt,delta)~X1+X2+X3, data=data, distribution=dist, weights=w, offset=rep(0, N),
                 train_params=params, var_monotone=c(0, 0, 0), keep_gbm_data=TRUE, cv_folds=5, is_verbose=FALSE)
    
    best_iter <- gbm_perf(gbm1, method="test") # returns test set estimate of best number of trees
  
    # make some new data
    set.seed(2)
    N <- 1000
    X1 <- runif(N)
    X2 <- runif(N)
    X3 <- factor(sample(letters[1:4],N,replace=T))
    mu <- c(-1,0,1,2)[as.numeric(X3)]

    f <- 0.5*sin(3*X1 + 5*X2^2 + mu/10)  # -0.5 <= f <= 0.5 via sin fn.
    tt.surv <- rexp(N,exp(f))
    tt.cens <- rexp(N,0.5)

    data2 <- data.frame(tt=apply(cbind(tt.surv,tt.cens),1,min),
                        delta=as.numeric(tt.surv <= tt.cens),
                        f=f,
                        X1=X1,X2=X2,X3=X3)

    # predict on the new data using "best" number of trees
    # f.predict will be on the canonical scale (logit,log,etc.)
    f.predict <- predict(gbm1, data2, best_iter)

    #plot(data2$f,f.predict)
    # Use observed sd
    expect_true(sd(data2$f - f.predict) < 0.4)
})
test_that("coxph - runs to completion with train_fraction of 1.0", {
  ## Needed packages
  require(survival)
  
  # Given data from the survival data
  # keep only certain baseline variables, subjects with longitudinal data
  temp <- subset(pbc, id%in%pbcseq$id, select=c(id:sex, stage))
  pbc1 <- tmerge(data1=temp, data2=temp, id=id, death = event(time, status))
  
  pbc2 <- tmerge(pbc1, pbcseq, id=id, ascites = tdc(day, ascites),
                 bili = tdc(day, bili), albumin = tdc(day, albumin),
                 protime = tdc(day, protime), alk.phos = tdc(day, alk.phos))
  # Training params
  params_1 <- training_params(interaction_depth = 3, num_train=round(1.0 * nrow(pbc)), num_trees = 500, id=seq(nrow(pbc)),
                              shrinkage=.01)
  params_2 <- training_params(interaction_depth = 3, num_train=round(1.0 * nrow(pbc2)), num_trees = 500, id=seq(nrow(pbc2)),
                              shrinkage=.01)
  
  # Then expect no errors when performing gbm fit with train_fraction = 1.0
  # GBM fit baseline data
  expect_error(gbmt(Surv(time, status==2) ~ bili + protime + albumin + alk.phos, data=pbc,
                    distribution=gbm_dist("CoxPH"), train_params=params_1), NA)
  
  # GBM fit using start/stop times to get time-dependent covariates
  expect_error(gbmt(Surv(tstart, tstop, death==2) ~ bili + protime + albumin + alk.phos, 
                   data=pbc2, distribution=gbm_dist("CoxPH"), train_params=params_2), NA)
})
test_that("coxph - runs to completion with train_fraction < 1.0 and cv_folds > 1", {
  ## Needed packages
  require(survival)
  
  # Given data from the survival data
  # keep only certain baseline variables, subjects with longitudinal data
  temp <- subset(pbc, id%in%pbcseq$id, select=c(id:sex, stage))
  pbc1 <- tmerge(data1=temp, data2=temp, id=id, death = event(time, status))
  
  pbc2 <- tmerge(pbc1, pbcseq, id=id, ascites = tdc(day, ascites),
                 bili = tdc(day, bili), albumin = tdc(day, albumin),
                 protime = tdc(day, protime), alk.phos = tdc(day, alk.phos))
  
  # Training params
  params_1 <- training_params(interaction_depth = 3, num_train=round(0.8 * nrow(pbc)), num_trees = 500, id=seq(nrow(pbc)),
                              shrinkage=.01)
  params_2 <- training_params(interaction_depth = 3, num_train=round(0.8 * nrow(pbc2)), num_trees = 500, id=seq(nrow(pbc2)),
                              shrinkage=.01)

  
  # Then expect no errors when performing gbm fit with train_fraction < 1.0 and cv_folds > 1
  # GBM fit baseline data
  expect_error(gbmt(Surv(time, status==2) ~ bili + protime + albumin + alk.phos, data=pbc,  distribution=gbm_dist("CoxPH"),
                   train_params=params_1), NA)
  expect_error(gbmt(Surv(time, status==2) ~ bili + protime + albumin + alk.phos, data=pbc,  distribution=gbm_dist("CoxPH"),
                   train_params=params_1, cv_folds=5), NA)
  
  # GBM fit using start/stop times to get time-dependent covariates
  expect_error(gbmt(Surv(tstart, tstop, death==2) ~ bili + protime + albumin + alk.phos, 
                   data=pbc2, distribution=gbm_dist("CoxPH"), train_params=params_2), NA)
  expect_error(gbmt(Surv(tstart, tstop, death==2) ~ bili + protime + albumin + alk.phos, 
                    data=pbc2, distribution=gbm_dist("CoxPH"), train_params=params_2, cv_folds=5), NA)
})
test_that("coxph cv_folds - runs to completion with start-stop, id'ed and stratified dataset", {
  ## Needed packages
  require(survival)
  
  # Given data from the survival package
  cgd2 <- cgd[cgd$enum==1,]
  
  # Training params
  params_1 <- training_params(interaction_depth = 1, num_train=nrow(cgd2), num_trees = 500, id=seq(nrow(cgd2)),
                            shrinkage=.01)
  params_2 <- training_params(interaction_depth = 1, num_train=round(0.8 * nrow(cgd2)), num_trees = 500, id=seq(nrow(cgd2)),
                              shrinkage=.01)
  params_3 <- training_params(interaction_depth = 3, num_train=round(1.0 * length(unique(cgd$id))), num_trees = 500, id=cgd$id,
                              shrinkage=.01)
  params_4 <- training_params(interaction_depth = 3, num_train=round(0.8 * length(unique(cgd$id))), num_trees = 500, id=cgd$id,
                              shrinkage=.01)
  
  # Then fitting a gbm model should throw no errors - with cv_folds > 1
  expect_error(gbmt(Surv(tstop, status) ~ age + sex + inherit +
                     steroids + propylac + hos.cat, data=cgd2, 
                   distribution = gbm_dist("CoxPH"), train_params=params_1, cv_folds=10), NA)
  expect_error(gbmt(Surv(tstop, status) ~ age + sex + inherit +
                     steroids + propylac + hos.cat, data=cgd2, 
                   distribution = gbm_dist("CoxPH"),  train_params=params_2, cv_folds=5), NA)
  
  expect_error(gbmt(Surv(tstart, tstop, status) ~ age + sex + inherit +
                     steroids + propylac, data=cgd,
                   distribution = gbm_dist("CoxPH", strata=cgd$hos.cat), train_params=params_3, cv_folds=10), NA)
  expect_error(gbmt(Surv(tstart, tstop, status) ~ age + sex + inherit +
                     steroids + propylac, data=cgd, 
                    distribution = gbm_dist("CoxPH", strata=cgd$hos.cat), train_params=params_4, cv_folds=10), NA)
})
test_that("bernoulli works", {
    set.seed(1)

    # create some data
    N <- 1000
    X1 <- runif(N)
    X2 <- runif(N)
    X3 <- factor(sample(letters[1:4],N,replace=T))
    mu <- c(-1,0,1,2)[as.numeric(X3)]

    p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
    Y <- rbinom(N,1,p)

    # random weights if you want to experiment with them
    w <- rexp(N)
    w <- N*w/sum(w)

    data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
    offset <- rep(0, N)
    
    # Set up for new API
    params <- training_params(num_trees=3000, interaction_depth=3, min_num_obs_in_node=10, 
                              shrinkage=0.001, bag_fraction=0.5, id=seq(nrow(data)), num_train=N/2, num_features=3)
    dist <- gbm_dist("Bernoulli")
    fit <- gbmt(Y~X1+X2+X3, data=data, distribution=dist, weights=w, offset=offset,
                train_params=params, var_monotone=c(0, 0, 0), keep_gbm_data=TRUE, cv_folds=5, is_verbose = FALSE)
    
    best_iter <- gbm_perf(fit, method="test") # returns test set estimate of best number of trees

    # make some new data
    set.seed(2)
    N <- 1000
    X1 <- runif(N)
    X2 <- runif(N)
    X3 <- factor(sample(letters[1:4],N,replace=T))
    mu <- c(-1,0,1,2)[as.numeric(X3)]

    p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
    Y <- rbinom(N,1,p)
    data2 <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)

    # predict on the new data using "best" number of trees
    # f.predict will be on the canonical scale (logit,log,etc.)
    f.1.predict <- predict(fit, data2, num_trees=best_iter)

    # compute quantity prior to transformation
    f.new = sin(3*X1) - 4*X2 + mu

    # Base the validation tests on observed discrepancies
    expect_true(sd(f.new - f.1.predict) < 1.0)
})
test_that("relative influence picks out true predictors", {
    set.seed(1234)
    X1 <- matrix(nrow=1000, ncol=50)
    X1 <- apply(X1, 2, function(x) rnorm(1000)) # Random noise
    X2 <- matrix(nrow=1000, ncol=5)
    X2 <- apply(X2, 2, function(x) c(rnorm(500), rnorm(500, 3))) # Real predictors
    cls <- rep(c(0, 1), ea=500) # Class
    X <- data.frame(cbind(X1, X2, cls))
    
    # Construct model
    params <- training_params(num_trees=1000, interaction_depth=2, min_num_obs_in_node=10, 
                              shrinkage=0.01, bag_fraction=0.5, id=seq(1000), num_train=1000, num_features=55)
    dist <- gbm_dist("Bernoulli")
    fit <- gbmt(cls~ ., data=X, distribution=dist, weights=rep(1, 1000), offset=rep(0, 1000),
                train_params=params, keep_gbm_data=TRUE, cv_folds=5, is_verbose = FALSE)
    
    # Check can get out real predictors
    ri <- relative_influence(fit, rescale=TRUE, sort_it=TRUE)
    wh <- names(ri)[1:5]
    res <- sum(wh %in% paste("V", 51:55, sep = ""))
    expect_equal(res, 5)
})
test_that("Conversion of 2 factor Y is successful", {
  
  NumY <- sample(c(0, 1), size=1000, replace=TRUE)
  FactY = factor(ifelse(NumY==1, "Yes", "No"), levels=c("No", "Yes"))

  PredX <-
  data.frame(
    x1 = runif(1000)
    ,x2 = runif(1000)
  )

  # Fit 1 - Rel Influence
  set.seed(32479)
  params <- training_params(num_trees=50, interaction_depth=3, min_num_obs_in_node=10, 
                            shrinkage=0.001, bag_fraction=0.5, id=seq(1000), num_train=1000, num_features=2)
  g1 <- gbmt(y ~ ., data = data.frame(y = NumY, PredX)
            , distribution = gbm_dist("Bernoulli"), train_params=params, is_verbose = FALSE)
  rig1 <- relative_influence(g1, num_trees=50)
  
  # Fit 2 - Rel Influence
   set.seed(32479)
  g2 <- gbmt(y ~ ., data = data.frame(y = FactY, PredX)
           ,  distribution = gbm_dist("Bernoulli"), train_params=params, is_verbose = FALSE)
   rig2 <- relative_influence(g2, num_trees=50)
  
  expect_equal(rig1, rig2)
})
test_that("Cross Validations group generation - Bernoulli", {
  
  NumY <- sample(rep(c(0,1), 500), size=1000)
  NumX <- matrix(rep(0, size=1000), nrow=1000, ncol=1)
  dist <- gbm_dist("Bernoulli")
  data <- gbm_data(NumX, NumY, weights=rep(1, length=1000), offset=rep(0, length=1000))
  
  cv_groups <- create_cv_groups(data, dist, training_params(num_train=1000, id=1:1000, num_features=1), 
                                cv_folds=4, cv_class_stratify=FALSE, fold_id=NULL)
  expect_true(all(table(cv_groups)==250))
  
  cv_stratified <- create_cv_groups(data, dist, training_params(num_train=1000, id=1:1000, num_features=1), 
                                    cv_folds=4, cv_class_stratify=TRUE, fold_id=NULL)
  
  Strats <- sapply(1:4, function(x){table(NumY[cv_stratified == 1])})
  
  expect_true(all(Strats == 125))
  
})

