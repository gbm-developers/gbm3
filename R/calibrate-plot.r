#' Calibration plot
#' 
#' An experimental diagnostic tool that plots the fitted values versus the
#' actual average values. Currently developed for only
#' \code{distribution="Bernoulli"}.
#' 
#' Uses natural splines to estimate E(y|p). Well-calibrated predictions imply
#' that E(y|p) = p. The plot also includes a pointwise 95% confidence band.
#' 
#' @param y the outcome 0-1 variable
#' @param p the predictions estimating E(y|x)
#' @param distribution the loss function used in creating \code{p}.
#' \code{bernoulli} and \code{poisson} are currently the only special options.
#' All others default to squared error assuming \code{gaussian}
#' @param replace determines whether this plot will replace or overlay the
#' current plot.  \code{replace=FALSE} is useful for comparing the calibration
#' of several methods
#' @param line.par graphics parameters for the line
#' @param shade.col color for shading the 2 SE region. \code{shade.col=NA}
#' implies no 2 SE region
#' @param shade.density the \code{density} parameter for \code{\link{polygon}}
#' @param rug.par graphics parameters passed to \code{\link{rug}}
#' @param xlab x-axis label corresponding to the predicted values
#' @param ylab y-axis label corresponding to the observed average
#' @param xlim,ylim x and y-axis limits. If not specified the function will
#' select limits
#' @param knots,df these parameters are passed directly to
#' \code{\link[splines]{ns}} for constructing a natural spline smoother for the
#' calibration curve
#' @param ...  other graphics parameters passed on to the plot function
#' @return \code{calibrate.plot} returns no values.
#' @author Greg Ridgeway \email{gregridgeway@@gmail.com}
#' @references J.F. Yates (1982). "External correspondence: decomposition of
#' the mean probability score," Organisational Behaviour and Human Performance
#' 30:132-156.
#' 
#' D.J. Spiegelhalter (1986). "Probabilistic Prediction in Patient Management
#' and Clinical Trials," Statistics in Medicine 5:421-433.
#' @keywords hplot
#' @examples
#' 
#' \dontrun{
#' library(rpart)
#' data(kyphosis)
#' y <- as.numeric(kyphosis$Kyphosis)-1
#' x <- kyphosis$Age
#' glm1 <- glm(y~poly(x,2),family=binomial)
#' p <- predict(glm1,type="response")
#' calibrate.plot(y, p, xlim=c(0,0.6), ylim=c(0,0.6))
#' }
#' @importFrom splines ns
calibrate.plot <- function(y, p,
                           distribution="Bernoulli",
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
      if(distribution=="gamma")
      {
         se.lower[se.lower < 0] <- 0
      }
      if(distribution=="tweedie")
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
