if(FALSE)
{
   library(gbm)

   N <- 10000
   X1 <- runif(N)
   X2 <- 2*runif(N)
   X3 <- ordered(sample(letters[1:4],N,replace=TRUE),levels=letters[4:1])
   X4 <- factor(sample(letters[1:6],N,replace=TRUE))
   X5 <- factor(sample(letters[1:3],N,replace=TRUE))
   X6 <- 3*runif(N)
   mu <- c(-1,0,1,2)[as.numeric(X3)]

   SNR <- 10 # signal-to-noise ratio
   Y <- X1**1.5 + 2 * (X2**.5) + mu
   sigma <- sqrt(var(Y)/SNR)
   Y <- Y + rnorm(N,0,sigma)

   # introduce some missing values
   X1[sample(1:N,size=500)] <- NA
   X4[sample(1:N,size=300)] <- NA

   data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)

   # fit initial model
   shrink <- c(0.1,0.05,0.01,0.005,0.001)
   err <- vector("list",length(shrink))

   for(i in 1:length(shrink))
   {
      gbm1 <- gbm(Y~X1+X2+X3+X4+X5+X6,
                     data=data,
                     distribution="gaussian",
                     n.trees=10000,
                     shrinkage=shrink[i],
                     interaction.depth=3,
                     bag.fraction = 0.5,
                     train.fraction = 0.2,
                     n.minobsinnode = 10,
                     keep.data=FALSE,
                     verbose=TRUE)
      err[[i]] <- gbm1$valid.error
   }

   ylim <- range(unlist(lapply(err,range)))
   mins <- min(sapply(err,min))
   ylim <- c(0.19,0.21)
   postscript("shrinkage-v-iterations.eps",horizontal=FALSE,width=9,height=6)
   plot(0,0,ylim=ylim,xlim=c(0,10000),type="n",xlab="Iterations",ylab="Squared error")
   for(i in 1:length(shrink))
   {
      x <- which.min(err[[i]])
      y <- err[[i]][x]
      j <- round(seq(1,10000,length=500))
      j <- sort(c(j,x))
      #k <- which((err[[i]][j] > ylim[1]) & (err[[i]][j] < ylim[2]))
      #k <- unique(c(1,k))
      k <- 1:length(j)
      lines(j[k],err[[i]][j][k],col=i)
      rug(x, col=i)
      text(x,y-0.0005,as.character(shrink[i]),adj=1)
   }
   abline(h=min(mins))
   dev.off()
}