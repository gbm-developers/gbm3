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
#' \code{Bernoulli} and \code{Poisson} are currently the only special options.
#' All others default to squared error assuming \code{Gaussian}
#' @param replace determines whether this plot will replace or overlay the
#' current plot.  \code{replace=FALSE} is useful for comparing the calibration
#' of several methods
#' @param line.par graphics parameters for the line
#' @param shade_col color for shading the 2 SE region. \code{shade_col=NA}
#' implies no 2 SE region
#' @param shade_density the \code{density} parameter for \code{\link{polygon}}
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
#' calibrate_plot(y, p, xlim=c(0,0.6), ylim=c(0,0.6))
#' }
#' @importFrom splines ns
calibrate_plot <- function(y, p,
                           distribution="Bernoulli",
                           replace=TRUE,
                           line.par=list(col="black"),
                           shade_col="lightyellow",
                           shade_density=NULL,
                           rug.par=list(side=1),
                           xlab="Predicted value",
                           ylab="Observed average",
                           xlim=NULL,ylim=NULL,
                           knots=NULL,df=6,
                           ...) {
  # Some initial checks
  if(is.null(knots) && is.null(df))
    stop("Either knots or df must be specified")
  if((df != round(df)) || (df<1))
    stop("df must be a positive integer")
  
  # Set up data and family
  data <- data.frame(y=y,p=p)
  family1 <- switch(tolower(distribution),
                    bernoulli=binomial,
                    poisson=poisson,
                    gaussian)
  
  # Fit a generalized linear model using splined predictions
  # as predictor variables
  gam1 <- glm(y~ns(p, df=df, knots=knots), data=data, family=family1)
  
  # Set up x and y values (new pre)
  x <- seq(min(p),max(p),length=200)
  glm_preds <- predict(gam1, newdata=data.frame(p=x), se.fit=TRUE, type="response")
  
  # Remove NAs
  x <- x[!is.na(glm_preds$fit)]
  glm_preds$se.fit <- glm_preds$se.fit[!is.na(glm_preds$fit)]
  glm_preds$fit <- glm_preds$fit[!is.na(glm_preds$fit)]
  
  # Set shading limits and 
  if(!is.na(shade_col)) {
    se.lower <- glm_preds$fit-2*glm_preds$se.fit
    se.upper <- glm_preds$fit+2*glm_preds$se.fit
    
    # Truncate upper and lower values - so they
    # make sense for the distribution under consideration
    if(tolower(distribution) %in% c("bernoulli", "poisson", "gamma", "tweedie"))
      se.lower[se.lower < 0] <- 0
    if(tolower(distribution)=="bernoulli") se.upper[se.upper > 1] <- 1
    
    # Set x and y limits
    if(is.null(xlim)) xlim <- range(se.lower,se.upper,x)
    if(is.null(ylim)) ylim <- range(se.lower,se.upper,x)
    
  } else {
    if(is.null(xlim)) xlim <- range(yy$fit,x)
    if(is.null(ylim)) ylim <- range(yy$fit,x)
  }
  
  
  # Plot
  if(replace) {
    plot(0,0,
         type="n",
         xlab=xlab,ylab=ylab,
         xlim=xlim,ylim=ylim,
         ...)
  }
  
  # Add shaded
  if(!is.na(shade_col)) {
    polygon(c(x,rev(x),x[1]),
            c(se.lower,rev(se.upper),se.lower[1]),
            col=shade_col,
            border=NA,
            density=shade_density)
  }
  lines(x, glm_preds$fit,col=line.par$col)
  quantile_rug(p,side=rug.par$side)
  abline(0, 1, col="red")
}
