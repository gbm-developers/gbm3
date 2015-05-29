#' Marginal plots of fitted gbm objects
#' 
#' Plots the marginal effect of the selected variables by "integrating" out the
#' other variables.
#' 
#' \code{plot.gbm} produces low dimensional projections of the
#' \code{\link{gbm.object}} by integrating out the variables not included in
#' the \code{i.var} argument. The function selects a grid of points and uses
#' the weighted tree traversal method described in Friedman (2001) to do the
#' integration. Based on the variable types included in the projection,
#' \code{plot.gbm} selects an appropriate display choosing amongst line plots,
#' contour plots, and \code{\link[lattice]{lattice}} plots. If the default
#' graphics are not sufficient the user may set \code{return.grid=TRUE}, store
#' the result of the function, and develop another graphic display more
#' appropriate to the particular example.
#' 
#' @param x a \code{\link{gbm.object}} fitted using a call to \code{\link{gbm}}
#' @param i.var a vector of indices or the names of the variables to plot. If
#' using indices, the variables are indexed in the same order that they appear
#' in the initial \code{gbm} formula.  If \code{length(i.var)} is between 1 and
#' 3 then \code{plot.gbm} produces the plots. Otherwise, \code{plot.gbm}
#' returns only the grid of evaluation points and their average predictions
#' @param n.trees the number of trees used to generate the plot. Only the first
#' \code{n.trees} trees will be used
#' @param continuous.resolution The number of equally space points at which to
#' evaluate continuous predictors
#' @param grid.levels A list containing the points at which to evaluate each
#' predictor in \code{i.var} (in the same order as \code{i.var}). For
#' continuous predictors this is usually a regular sequence of values within
#' the range of the variable. For categorical predictors, the points are the
#' levels of the factor. When \code{length(i.var)} is one, the values can be
#' provided directly, outside a list.  This is NULL by default and generated
#' automatically from the data, using \code{continuous.resolution} for
#' continuous predictors. Forcing the values can be useful to evaluate two
#' models on the same exact range
#' @param return.grid if \code{TRUE} then \code{plot.gbm} produces no graphics
#' and only returns the grid of evaluation points and their average
#' predictions. This is useful for customizing the graphics for special
#' variable types or for dimensions greater than 3
#' @param type the type of prediction to plot on the vertical axis. See
#' \code{predict.gbm}
#' @param \dots other arguments passed to the plot function
#' @return Nothing unless \code{return.grid} is true then \code{plot.gbm}
#' produces no graphics and only returns the grid of evaluation points and
#' their average predictions.
#' @author Greg Ridgeway \email{gregridgeway@@gmail.com}
#' @seealso \code{\link{gbm}}, \code{\link{gbm.object}},
#' \code{\link[graphics]{plot}}
#' @references J.H. Friedman (2001). "Greedy Function Approximation: A Gradient
#' Boosting Machine," Annals of Statistics 29(4).
#' @keywords hplot
#' @export
plot.gbm <- function(x,
                     i.var=1,
                     n.trees=x$n.trees,
                     continuous.resolution=100,
                     grid.levels=NULL,
                     return.grid=FALSE,
                     type="link",
                     ...)
{
   if (!is.element(type, c("link", "response"))){
      stop( "type must be either 'link' or 'response'")
   }

   if(all(is.character(i.var)))
   {
      i <- match(i.var,x$var.names)
      if(any(is.na(i)))
      {
         stop("Plot variables not used in gbm model fit: ",i.var[is.na(i)])
      } else
      {
         i.var <- i
      }
   }

   if((min(i.var)<1) || (max(i.var)>length(x$var.names)))
   {
      warning("i.var must be between 1 and ",length(x$var.names))
   }
   if(n.trees > x$n.trees)
   {
      warning(paste("n.trees exceeds the number of trees in the model, ",x$n.trees,
                    ". Plotting using ",x$n.trees," trees.",sep=""))
      n.trees <- x$n.trees
   }

   if(length(i.var) > 3)
   {
      warning("gbm.int.plot creates up to 3-way interaction plots.\nplot.gbm will only return the plotting data structure.")
      return.grid = TRUE
   }

   # generate grid to evaluate gbm model
   if (is.null(grid.levels)) {
      grid.levels <- vector("list",length(i.var))
      for(i in 1:length(i.var))
      {
        # continuous
        if(is.numeric(x$var.levels[[i.var[i]]]))
        {
           grid.levels[[i]] <- seq(min(x$var.levels[[i.var[i]]]),
                                   max(x$var.levels[[i.var[i]]]),
                                   length=continuous.resolution)
        }
        # categorical or ordered
        else
        {
           grid.levels[[i]] <- as.numeric(factor(x$var.levels[[i.var[i]]],
                                                 levels=x$var.levels[[i.var[i]]]))-1
        }
      }
   }
   else
   {
      # allow grid.levels to not be a list when there is only one predictor
      if (length(i.var) == 1 & !is.list(grid.levels)) {
         grid.levels <- list(grid.levels)
      }
      # check compatibility in length, at least
      if (length(grid.levels) != length(i.var)) {
         stop("Need grid.levels for all variables in i.var")
      }
      # convert levels for categorical predictors into numbers
      for(i in 1:length(i.var)) {
         if(!is.numeric(x$var.levels[[i.var[i]]]))
         {
            grid.levels[[i]] <- as.numeric(grid.levels[[i]]) - 1
         }
      }
   }

   X <- expand.grid(grid.levels)
   names(X) <- paste("X",1:length(i.var),sep="")

   # Next if block for compatibility with objects created with 1.6
   if (is.null(x$num.classes)){
       x$num.classes <- 1
   }

   # evaluate at each data point
   y <- .Call("gbm_plot",
              X = data.matrix(X),
              i.var = as.integer(i.var-1),
              n.trees = as.integer(n.trees) ,
              initF = as.double(x$initF),
              trees = x$trees,
              c.splits = x$c.splits,
              var.type = as.integer(x$var.type),
              PACKAGE = "gbm")

   if(is.element(x$distribution$name, c("bernoulli", "pairwise")) && type=="response") {
      X$y <- 1/(1+exp(-y))
   }
   else if ((x$distribution$name=="poisson") && (type=="response")){
      X$y <- exp(y)
   }
   else if ((x$distribution$name=="gamma") && (type=="response")){
      X$y <- exp(y)
   }
   else if ((x$distribution$name=="tweedie") && (type=="response")){
      X$y <- exp(y)
   }
   else if (type=="response"){
      warning("type 'response' only implemented for 'bernoulli', 'poisson', 'gamma', 'tweedie', and 'pairwise'. Ignoring" )
   }
   else { X$y <- y }

   # transform categorical variables back to factors
   f.factor <- rep(FALSE,length(i.var))
   for(i in 1:length(i.var))
   {
      if(!is.numeric(x$var.levels[[i.var[i]]]))
      {
         X[,i] <- factor(x$var.levels[[i.var[i]]][X[,i]+1],
                         levels=x$var.levels[[i.var[i]]])
         f.factor[i] <- TRUE
      }
   }

   if(return.grid)
   {
      names(X)[1:length(i.var)] <- x$var.names[i.var]
      return(X)
   }

   # create the plots
   if(length(i.var)==1)
   {
      if(!f.factor)
      {
         j <- order(X$X1)

         if (is.element(x$distribution$name, c("bernoulli", "pairwise"))) {
            if ( type == "response" ){
               ylabel <- "Predicted probability"
            }
            else {
               ylabel <- paste("f(",x$var.names[i.var],")",sep="")
            }
            plot( X$X1, X$y , type = "l", xlab = x$var.names[i.var], ylab=ylabel )
         }
         else if ( x$distribution$name == "poisson" ){
            if (type == "response" ){
               ylabel <- "Predicted count"
            }
            else{
               ylabel <- paste("f(",x$var.names[i.var],")",sep="")
            }
            plot( X$X1, X$y , type = "l", xlab = x$var.names[i.var], ylab=ylabel )
         }
         else {
            plot(X$X1,X$y,
                 type="l",
                 xlab=x$var.names[i.var],
                 ylab=paste("f(",x$var.names[i.var],")",sep=""),...)
         }
      }
      else
      {
         if (is.element(x$distribution$name, c("bernoulli", "pairwise")) && type == "response" ){
            ylabel <- "Predicted probability"
            plot( X$X1, X$y, type = "l", xlab=x$var.names[i.var], ylab=ylabel )
         }
         else if ( x$distribution$name == "poisson" & type == "response" ){
            ylabel <- "Predicted count"
            plot( X$X1, X$y, type = "l", xlab=x$var.names[i.var], ylab=ylabel )
         }
         else {
            plot(X$X1,X$y,
                 type="l",
                 xlab=x$var.names[i.var],
                 ylab=paste("f(",x$var.names[i.var],")",sep=""),...)
         }
      }
   }
   else if(length(i.var)==2)
   {
      if(!f.factor[1] && !f.factor[2])
      {
         print(levelplot(y~X1*X2,data=X,
                         xlab=x$var.names[i.var[1]],
                         ylab=x$var.names[i.var[2]],...))
         
      }
      else if(f.factor[1] && !f.factor[2])
      {
         print(xyplot(y~X2|X1,data=X,
                      xlab=x$var.names[i.var[2]],
                      ylab=paste("f(",x$var.names[i.var[1]],",",x$var.names[i.var[2]],")",sep=""),
                      type="l",
                      panel = panel.xyplot,
                      ...))
      }
      else if(!f.factor[1] && f.factor[2])
      {
         print(xyplot(y~X1|X2,data=X,
                      xlab=x$var.names[i.var[1]],
                      ylab=paste("f(",x$var.names[i.var[1]],",",x$var.names[i.var[2]],")",sep=""),
                      type="l",
                      panel = panel.xyplot,
                      ...))
     }
      else
      {
         print(stripplot(X1~y|X2,data=X,
                         xlab=x$var.names[i.var[2]],
                         ylab=paste("f(",x$var.names[i.var[1]],",",x$var.names[i.var[2]],")",sep=""),
                         ...))
         
     }
   }
   else if(length(i.var)==3)
   {
      i <- order(f.factor)
      X.new <- X[,i]
      X.new$y <- X$y
      names(X.new) <- names(X)

      # 0 factor, 3 continuous
      if(sum(f.factor)==0)
      {
         X.new$X3 <- equal.count(X.new$X3)
         print(levelplot(y~X1*X2|X3,data=X.new,
                         xlab=x$var.names[i.var[i[1]]],
                         ylab=x$var.names[i.var[i[2]]],...))
         
      }
      # 1 factor, 2 continuous
      else if(sum(f.factor)==1)
      {
         print(levelplot(y~X1*X2|X3,data=X.new,
                         xlab=x$var.names[i.var[i[1]]],
                         ylab=x$var.names[i.var[i[2]]],...))
         
      }
      # 2 factors, 1 continuous
      else if(sum(f.factor)==2)
      {
         print(xyplot(y~X1|X2*X3,data=X.new,
                      type="l",
                      xlab=x$var.names[i.var[i[1]]],
                      ylab=paste("f(",paste(x$var.names[i.var[1:3]],collapse=","),")",sep=""),
                      panel = panel.xyplot,
                      ...))
      }
      # 3 factors, 0 continuous
      else if(sum(f.factor)==3)
      {
         print(stripplot(X1~y|X2*X3,data=X.new,
                         xlab=x$var.names[i.var[i[1]]],
                         ylab=paste("f(",paste(x$var.names[i.var[1:3]],collapse=","),")",sep=""),
                         ...))

      }
   }
}
