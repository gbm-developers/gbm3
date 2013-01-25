quantile.rug <- function(x,prob=(0:10)/10,...)
{
   quants <- quantile(x[!is.na(x)],prob=prob)
   if(length(unique(quants)) < length(prob))
   {
      quants <- jitter(quants)
   }
   rug(quants,...)
}

calibrate.plot <- function(y,p,
                           distribution="bernoulli",
                           replace=TRUE,
                           line.par=list(col="black"),
                           shade.col="lightyellow",
                           shade.density=NULL,
                           rug.par=list(side=1),
                           xlab="Predicted value",
                           ylab="Observed average",
                           xlim=NULL,ylim=NULL,
                           knots=NULL,df=6,
                           ...)
{
   data <- data.frame(y=y,p=p)

   if(is.null(knots) && is.null(df))
      stop("Either knots or df must be specified")
   if((df != round(df)) || (df<1))
      stop("df must be a positive integer")

   if(distribution=="bernoulli")
   {
      family1 = binomial
   } else if(distribution=="poisson")
   {
      family1 = poisson
   } else
   {
      family1 = gaussian
   }
   gam1 <- glm(y~ns(p,df=df,knots=knots),data=data,family=family1)

   x <- seq(min(p),max(p),length=200)
   yy <- predict(gam1,newdata=data.frame(p=x),se.fit=TRUE,type="response")

   x <- x[!is.na(yy$fit)]
   yy$se.fit <- yy$se.fit[!is.na(yy$fit)]
   yy$fit <- yy$fit[!is.na(yy$fit)]

   if(!is.na(shade.col))
   {
      se.lower <- yy$fit-2*yy$se.fit
      se.upper <- yy$fit+2*yy$se.fit
      if(distribution=="bernoulli")
      {
         se.lower[se.lower < 0] <- 0
         se.upper[se.upper > 1] <- 1
      }
      if(distribution=="poisson")
      {
         se.lower[se.lower < 0] <- 0
      }
      if(is.null(xlim)) xlim <- range(se.lower,se.upper,x)
      if(is.null(ylim)) ylim <- range(se.lower,se.upper,x)
   }
   else
   {
      if(is.null(xlim)) xlim <- range(yy$fit,x)
      if(is.null(ylim)) ylim <- range(yy$fit,x)
   }
   if(replace)
   {
      plot(0,0,
           type="n",
           xlab=xlab,ylab=ylab,
           xlim=xlim,ylim=ylim,
           ...)
   }
   if(!is.na(shade.col))
   {
      polygon(c(x,rev(x),x[1]),
              c(se.lower,rev(se.upper),se.lower[1]),
              col=shade.col,
              border=NA,
              density=shade.density)
   }
   lines(x,yy$fit,col=line.par$col)
   quantile.rug(p,side=rug.par$side)
   abline(0,1,col="red")
}
