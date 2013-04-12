set.seed(06182001)

# number of replicates
n.reps <- 20
# should data be loaded from the web? If FALSE use alt.path
load.from.web <- TRUE
run.all <- TRUE
# if data not downloaded from the web, give path to datasets
alt.path <- ""


n.datasets <- 12 # needs to match the number of datasets
i.data <- 0

squared.error.loss <- function(y,f.x)
{
    mean((y-f.x)^2)
}

bernoulli.loglikelihood <- function(y,f.x)
{
    mean(y*f.x - log(1+exp(f.x)))
}


if(run.all)
{
   dataset <- vector("list",n.datasets)

   # abalone
   i.data <- i.data + 1
   dataset[[i.data]] <-
      list(name="Abalone",
            distribution="gaussian",
            urlpath="http://ftp.ics.uci.edu/pub/machine-learning-databases/abalone/",
            filename="abalone.data",
            var.names=c("sex","length","diameter","height","whole.weight",
                        "shucked.weight","viscera.weight","shell.weight",
                        "Rings"),
            outcome="Rings",
            factors="sex",
            na.strings="",
            sep=",",
            shrinkage=0.02)

   # Adult
   i.data <- i.data + 1
   dataset[[i.data]] <-
      list(name="Adult",
            distribution="bernoulli",
            urlpath="http://ftp.ics.uci.edu/pub/machine-learning-databases/adult/",
            filename="adult.data",
            var.names=c("age","workclass","w","education","education.num",
                        "marital.status","occupation","relationship","race",
                        "male","capital.gain","capital.loss",
                        "hours.per.week","native.country","income"),
            outcome="income",
            factors=c("workclass","education","marital.status","occupation",
                     "relationship","race","native.country","male"),
            na.strings="?",
            sep=",",
            shrinkage=0.04)

   # Housing
   i.data <- i.data + 1
   dataset[[i.data]] <-
      list(name="Boston housing",
            distribution="gaussian",
            urlpath="http://ftp.ics.uci.edu/pub/machine-learning-databases/housing/",
            filename="housing.data",
            var.names=c("CRIM","ZN","INDUS","CHAS","NOX","RM","AGE",
                        "DIS","RAD","TAX","PTRATIO","B","LSTAT","MEDV"),
            factors=NULL,
            outcome="MEDV",
            na.strings="",
            sep="",
            shrinkage=0.005)

   # mushrooms
   i.data <- i.data + 1
   dataset[[i.data]] <-
      list(name="Mushrooms",
            distribution="bernoulli",
            urlpath="http://ftp.ics.uci.edu/pub/machine-learning-databases/mushroom/",
            filename="agaricus-lepiota.data",
            var.names=c("poisonous","cap-shape","cap-surface","cap-color",
                        "bruises","odor","gill-attachment",
                        "gill-spacing","gill-size","gill-color",
                        "stalk-shape","stalk-root","stalk-surface-above-ring",
                        "stalk-surface-below-ring","stalk-color-above-ring",
                        "stalk-color-below-ring","veil-type","veil-color",
                        "ring-number","ring-type","spore-print-color",
                        "population","habitat"),
            factors=c("cap-shape","cap-surface","cap-color",
                     "bruises","odor","gill-attachment",
                     "gill-spacing","gill-size","gill-color",
                     "stalk-shape","stalk-root","stalk-surface-above-ring",
                     "stalk-surface-below-ring","stalk-color-above-ring",
                     "stalk-color-below-ring","veil-type","veil-color",
                     "ring-number","ring-type","spore-print-color",
                     "population","habitat"),
            outcome="poisonous",
            drop.vars=c("veil-type"),
            na.strings="?",
            sep=",",
            shrinkage=0.05)

   # autoprices 1
   i.data <- i.data + 1
   dataset[[i.data]] <-
      list(name="Auto Prices",
            distribution="gaussian",
            urlpath="http://ftp.ics.uci.edu/pub/machine-learning-databases/autos/",
            filename="imports-85.data",
            var.names=c("symboling","normalizedlosses","make","fueltype",
                        "aspiration","ndoors","bodystyle",
                        "drivewheels","enginelocation", "wheelbase", "length",
                        "width", "height", "curbweight", "enginetype",
                        "numerofcylinders", "enginesize", "fuelsystem", "bore",
                        "stroke", "compressionratio", "horsepower", "peakrpm",
                        "citympg", "highwatmpg", "price"),
            factors=c("symboling","make","fueltype","aspiration","ndoors",
                     "bodystyle","drivewheels","enginelocation", "enginetype",
                     "numerofcylinders", "fuelsystem"),
            outcome="price",
            na.strings="?",
            sep=",",
            shrinkage=0.002)

   # auto MPG
   i.data <- i.data + 1
   dataset[[i.data]] <-
      list(name="Auto MPG",
            distribution="gaussian",
            urlpath="http://ftp.ics.uci.edu/pub/machine-learning-databases/auto-mpg/",
            filename="auto-mpg.data",
            var.names=c("mpg","cylinders","displacement","horsepower","weight",
                        "acceleration","modelyear","origin","carname"),
            factors=c("cylinders", "modelyear", "origin"),
            outcome="mpg",
            drop.vars=c("carname"),
            na.strings="?",
            sep="",
            shrinkage=0.005)

   # CPU
   i.data <- i.data + 1
   dataset[[i.data]] <-
      list(name="CPU Performance",
            distribution="gaussian",
            urlpath="http://ftp.ics.uci.edu/pub/machine-learning-databases/cpu-performance/",
            filename="machine.data",
            var.names=c("vendorname","modelname","myct","mmin","mmax",
                        "cach","chmin","chmax","prp","ERP"),
            factors=c("vendorname","modelname"),
            outcome="prp",
            na.strings="",
            drop.vars=c("vendorname","modelname"),
            sep=",",
            shrinkage=0.01)

   # credit
   i.data <- i.data + 1
   dataset[[i.data]] <-
      list(name="Credit rating",
            distribution="bernoulli",
            urlpath="http://ftp.ics.uci.edu/pub/machine-learning-databases/credit-screening/",
            filename="crx.data",
            var.names=c("A1","A2","A3","A4","A5","A6","A7","A8","A9","A10","A11",
                        "A12", "A13", "A14", "A15","CLASS"),
            factors=c("A1","A4", "A5", "A6", "A7", "A9", "A10", "A12", "A13","CLASS"),
            outcome="CLASS",
            na.strings="?",
            sep=",",
            shrinkage=0.005)

   # Haberman
   i.data <- i.data + 1
   dataset[[i.data]] <-
      list(name="Haberman",
            distribution="bernoulli",
            urlpath="http://ftp.ics.uci.edu/pub/machine-learning-databases/haberman/",
            filename="haberman.data",
            var.names=c("age","year","nodes","CLASS"),
            outcome="CLASS",
            factors=c("CLASS"),
            na.strings="",
            sep=",",
            shrinkage=0.001)

   # Diabetes
   i.data <- i.data + 1
   dataset[[i.data]] <-
      list(name="Diabetes",
            distribution="bernoulli",
            urlpath="http://ftp.ics.uci.edu/pub/machine-learning-databases/pima-indians-diabetes/",
            filename="pima-indians-diabetes.data",
            var.names=c("n_preg","plasma","blood-pre","triceps","serum",
                     "mass-index","pedigree","age","CLASS"),
            factors=c("CLASS"),
            outcome="CLASS",
            na.strings="?",
            sep=",",
            shrinkage=0.005)

   # Ionosphere
   i.data <- i.data + 1
   dataset[[i.data]] <-
      list(name="Ionosphere",
            distribution="bernoulli",
            urlpath="http://ftp.ics.uci.edu/pub/machine-learning-databases/ionosphere/",
            filename="ionosphere.data",
            var.names=c("A1","A2","A3","A4","A5","A6","A7","A8","A9","A10","A11",
                        "A12","A13","A14","A15","A16","A17","A18","A19","A20",
                        "A21","A22","A23","A24","A25","A26","A27","A28","A29",
                        "A30","A31","A32","A33","A34","CLASS"),
            factors=c("CLASS"),
            outcome="CLASS",
            na.strings="",
            sep=",",
            shrinkage=0.005)

   # Breast cancer
   i.data <- i.data + 1
   dataset[[i.data]] <-
      list(name="breast cancer",
            distribution="bernoulli",
            urlpath="http://ftp.ics.uci.edu/pub/machine-learning-databases/breast-cancer-wisconsin/",
            filename="breast-cancer-wisconsin.data",
            var.names=c("CODE","thickness","cellsize","cellshape","adhension",
                        "singleecell","bnuclei","chromatin","nnucleo","mitoses",
                        "CLASS"),
            factors=c("CODE","CLASS"),
            outcome="CLASS",
            drop.vars=c("CODE"),
            na.strings="?",
            sep=",",
            shrinkage=0.005)

   if(FALSE) # this dataset is not public, can substitute other datasets
   {
      # time in treatment
      i.data <- i.data + 1
      dataset[[i.data]] <-
         list(name="time in treatment",
               distribution="gaussian",
               urlpath="./",
               filename="txdet.csv",
               var.names=NULL,
               factors=c("b1","xsite4","b3new","b8new","s1a1new","m3dnew","e1new","e13anew"),
               outcome="txdet",
               drop.vars=c("xpid","xobs","maxcefu","recovfu","nontxdet","s7e5","r2f",
                           "r3a9","e4a6","l5p","v2.4","v2.7","v2.8"),
               na.strings="NA",
               sep=",",
               shrinkage=0.0022)
   }

   # Load datasets
   for(i.data in 1:n.datasets)
   # for(i.data in which(sapply(dataset,function(x){is.null(x$oob.iter)})))
   {
      # Progress
      cat("Dataset ",i.data,":",dataset[[i.data]]$name," N = ")
      filename <- paste(switch(load.from.web+1,
                              alt.path,
                              dataset[[i.data]]$url),
                        dataset[[i.data]]$filename,
                        sep="")
      dataset[[i.data]]$data <-
         read.table(file=filename,
                     na.strings=dataset[[i.data]]$na.strings,
                     sep=dataset[[i.data]]$sep,
                     header=is.null(dataset[[i.data]]$var.names))
      if(!is.null(dataset[[i.data]]$var.names))
      {
         names(dataset[[i.data]]$data) <- dataset[[i.data]]$var.names
      }

      # take care of nominal predictors
      for(j in dataset[[i.data]]$factors)
      {
         dataset[[i.data]]$data[,j] <- factor(dataset[[i.data]]$data[,j])
      }

      # take care of factor binary outcomes
      if( with(dataset[[i.data]],
         (distribution=="bernoulli") && is.factor(data[,outcome])) )
      {
         dataset[[i.data]]$data[,dataset[[i.data]]$outcome] <-
               with(dataset[[i.data]], as.numeric(data[,outcome])-1)
      }

      # drop observations with missing outcomes
      i <- with(dataset[[i.data]], !is.na(data[,outcome]))
      dataset[[i.data]]$data <- dataset[[i.data]]$data[i,]

      # drop selected predictor variables
      if(!is.null(dataset[[i.data]]$drop.vars))
      {
         j <- match(dataset[[i.data]]$drop.vars,names(dataset[[i.data]]$data))
         dataset[[i.data]]$data <- dataset[[i.data]]$data[,-j]
      }

      dataset[[i.data]]$loss <- switch(dataset[[i.data]]$distribution,
                                       gaussian=squared.error.loss,
                                       bernoulli=bernoulli.loglikelihood)

      cat(nrow(dataset[[i.data]]$data),"\n")
   }

   save(dataset,file="dataset.RData")
} # run.all


# make sure gbm is installed
if(!is.element("gbm",installed.packages()[,1]))
{
    stop("The gbm package is not installed. Use install.packages(\"gbm\") to install. On Unix machines this must be executed in an R session started as root or installed to a local library, see help(install.packages)")
}

library(gbm)
# loop over all the datasets
i.datasets <- which(sapply(dataset,function(x){is.null(x$oob.loss)}))
for(i.data in i.datasets)
# for(i.data in which(sapply(dataset,function(x){is.null(x$oob.iter)})))
{
   N <- nrow(dataset[[i.data]]$data)
   # Progress
   cat("Dataset ",i.data,":",dataset[[i.data]]$name," N = ",N,"\n",sep="")

   # construct model formula for this dataset
   formula.fit <- formula(paste(dataset[[i.data]]$outcome,"~ ."))

   # initialize prediction
   pred.oob <- pred.base <- pred.test33 <- pred.test20 <- pred.cv5 <- rep(0,N)

   # track iteration estimates
   dataset[[i.data]]$oob.iter    <- rep(NA,n.reps)
   dataset[[i.data]]$test33.iter <- rep(NA,n.reps)
   dataset[[i.data]]$test20.iter <- rep(NA,n.reps)
   dataset[[i.data]]$cv5.iter    <- rep(NA,n.reps)

   # do replicates
   for(i.rep in 1:n.reps)
   {
      cat("rep:",i.rep,"")
      i.train <- sample(1:N,size=0.75*N,replace=FALSE)
      i.valid <- (1:N)[-i.train]

      # use out-of-bag method
      cat("OOB, ")
      gbm1 <- gbm(formula.fit,
                  data=dataset[[i.data]]$data[i.train,],
                  distribution=dataset[[i.data]]$distribution,
                  train.fraction=1.0,
                  bag.fraction=0.5,
                  shrinkage=dataset[[i.data]]$shrinkage,
                  n.trees=1000,
                  verbose = FALSE)

      best.iter.oob <- gbm.perf(gbm1,method="OOB",plot.it=FALSE)
      while((gbm1$n.trees-best.iter.oob < 1000) &&
            !all(gbm1$oobag.improve[(gbm1$n.trees-100):gbm1$n.trees] < 1e-6))
      {
         gbm1 <- gbm.more(gbm1,1000)
         best.iter.oob <- gbm.perf(gbm1,method="OOB",plot.it=FALSE)
      }
      pred.oob[i.valid] <- predict(gbm1,
                                    newdata=dataset[[i.data]]$data[i.valid,],
                                    n.trees=best.iter.oob)
      dataset[[i.data]]$oob.iter[i.rep] <- best.iter.oob

      # use a 1/3 test set
      cat("33% test data, ")
      gbm1 <- gbm(formula.fit,
                  data=dataset[[i.data]]$data[i.train,],
                  distribution=dataset[[i.data]]$distribution,
                  train.fraction=2/3,
                  bag.fraction=0.5,
                  shrinkage=dataset[[i.data]]$shrinkage,
                  n.trees=1000,
                  verbose = FALSE)
      best.iter.test <- gbm.perf(gbm1,method="test",plot.it=FALSE)
      while((gbm1$n.trees-best.iter.test < 1000) &&
            !all(abs(gbm1$valid.error[(gbm1$n.trees-100):gbm1$n.trees]) < 1e-6))
      {
         gbm1 <- gbm.more(gbm1,1000)
         best.iter.test <- gbm.perf(gbm1,method="test",plot.it=FALSE)
      }
      pred.test33[i.valid] <- predict(gbm1,
                                    newdata=dataset[[i.data]]$data[i.valid,],
                                    n.trees=best.iter.test)
      dataset[[i.data]]$test33.iter[i.rep] <- best.iter.test

      # use a 20% test set
      cat("20% test data, ")
      gbm1 <- gbm(formula.fit,
                  data=dataset[[i.data]]$data[i.train,],
                  distribution=dataset[[i.data]]$distribution,
                  train.fraction=0.80,
                  bag.fraction=0.5,
                  shrinkage=dataset[[i.data]]$shrinkage,
                  n.trees=1000,
                  verbose = FALSE)
      best.iter.test <- gbm.perf(gbm1,method="test",plot.it=FALSE)
      while((gbm1$n.trees-best.iter.test < 1000) &&
            !all(abs(gbm1$valid.error[(gbm1$n.trees-100):gbm1$n.trees]) < 1e-6))
      {
         gbm1 <- gbm.more(gbm1,1000)
         best.iter.test <- gbm.perf(gbm1,method="test",plot.it=FALSE)
      }
      pred.test20[i.valid] <-
         predict(gbm1,
                  newdata=dataset[[i.data]]$data[i.valid,],
                  n.trees=best.iter.test)
      dataset[[i.data]]$test20.iter[i.rep] <- best.iter.test

      # use 5-fold cross-validation
      cat("5-fold CV")
      n.cv <- 5
      cv.group  <- sample(rep(1:n.cv,length=length(i.train)))
      max.iters <- round(best.iter.test*1.2)
      cv.loss   <- matrix(0,ncol=n.cv,nrow=max.iters)
      for(i.cv in 1:n.cv)
      {
         cat(".")
         i <- order(cv.group==i.cv) # used to put the held out obs last
         gbm1 <- gbm(formula.fit,
                     data=dataset[[i.data]]$data[i.train[i],],
                     distribution=dataset[[i.data]]$distribution,
                     train.fraction=mean(cv.group!=i.cv),
                     bag.fraction=0.5,
                     shrinkage=dataset[[i.data]]$shrinkage,
                     n.trees=max.iters,
                     verbose = FALSE)
         cv.loss[,i.cv] <- gbm1$valid.error
      }
      cat("\n")

      best.iter.cv <- which.min(apply(cv.loss,1,weighted.mean,w=table(cv.group)))
      gbm1 <- gbm(formula.fit,
                  data=dataset[[i.data]]$data[i.train,],
                  distribution=dataset[[i.data]]$distribution,
                  train.fraction=1.0,
                  bag.fraction=0.5,
                  shrinkage=dataset[[i.data]]$shrinkage,
                  n.trees=best.iter.cv,
                  verbose = FALSE)
      pred.cv5[i.valid] <-
         predict(gbm1,
                 newdata=dataset[[i.data]]$data[i.valid,],
                 n.trees=best.iter.cv)
      dataset[[i.data]]$cv5.iter[i.rep] <- best.iter.cv

      # baseline prediction
      pred.base[i.valid] <- gbm1$initF

      # evalute the methods
      dataset[[i.data]]$base.loss[i.rep] <-
         with(dataset[[i.data]], loss(data[i.valid,outcome],pred.base[i.valid]))
      dataset[[i.data]]$oob.loss[i.rep]  <-
         with(dataset[[i.data]], loss(data[i.valid,outcome],pred.oob[i.valid]))
      dataset[[i.data]]$test33.loss[i.rep] <-
         with(dataset[[i.data]], loss(data[i.valid,outcome],pred.test33[i.valid]))
      dataset[[i.data]]$test20.loss[i.rep] <-
         with(dataset[[i.data]], loss(data[i.valid,outcome],pred.test20[i.valid]))
      dataset[[i.data]]$cv5.loss[i.rep] <-
         with(dataset[[i.data]], loss(data[i.valid,outcome],pred.cv5[i.valid]))

      with(dataset[[i.data]],
          cat(oob.iter[i.rep],test33.iter[i.rep],test20.iter[i.rep],
              cv5.iter[i.rep],"\n"))
   }

   save.image(compress=TRUE)
}

#rm(dataset)
save.image(compress=TRUE)

results <- data.frame(problem=sapply(dataset,function(x){x$name}),
                      N=sapply(dataset,function(x){nrow(x$data)}),
                      d=sapply(dataset,function(x){ncol(x$data)-1}),
                      loss=sapply(dataset,function(x){x$distribution}),
                      base=sapply(dataset,function(x){mean(x$base.loss)}),
                      oob=sapply(dataset,function(x){mean(x$oob.loss)}),
                      test33=sapply(dataset,function(x){mean(x$test33.loss)}),
                      test20=sapply(dataset,function(x){mean(x$test20.loss)}),
                      cv5=sapply(dataset,function(x){mean(x$cv5.loss)}))
j <- match(c("base","oob","test33","test20","cv5"),names(results))
results[results$loss=="bernoulli",j] <- -2*results[results$loss=="bernoulli",j]
results$win <- c("base","oob","test33","test20","cv5")[apply(results[,j],1,which.min)]
results$oob.rank <- apply(results[,j],1,rank)[2,]
results$perf <- (results$base-results$oob)/apply(results$base-results[,j],1,max)

plot(0,0,ylim=c(0,14000),xlim=c(0,n.datasets+1),
     xlab="Dataset",ylab="Number of iterations",
     type="n",axes=FALSE)
lines(sapply(dataset,function(x){mean(x$oob.iter)}),    col="blue")
lines(sapply(dataset,function(x){mean(x$test33.iter)}), col="red")
lines(sapply(dataset,function(x){mean(x$test20.iter)}), col="green")
lines(sapply(dataset,function(x){mean(x$cv5.iter)}),    col="purple")
axis(1,at=1:n.datasets,labels=as.character(results$problem))
