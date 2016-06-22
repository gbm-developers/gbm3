#' Marginal plots of fitted gbm objects
#' 
#' Plots the marginal effect of the selected variables by "integrating" out the
#' other variables.
#' 
#' \code{plot.GBMFit} produces low dimensional projections of the
#' \code{\link{GBMFit}} object, see \code{\link{gbm2}}, by integrating out the variables not included in
#' the \code{var_index} argument. The function selects a grid of points and uses
#' the weighted tree traversal method described in Friedman (2001) to do the
#' integration. Based on the variable types included in the projection,
#' \code{plot.GBMFit} selects an appropriate display choosing amongst line plots,
#' contour plots, and \code{\link[lattice]{lattice}} plots. If the default
#' graphics are not sufficient the user may set \code{return.grid=TRUE}, store
#' the result of the function, and develop another graphic display more
#' appropriate to the particular example.
#' 
#' @param gbm_fit_obj a \code{GBMFit} object fitted using a call to \code{\link{gbm2}}
#' @param var_index a vector of indices or the names of the variables to plot. If
#' using indices, the variables are indexed in the same order that they appear
#' in the initial \code{gbm2} formula.  If \code{length(var_index)} is between 1 and
#' 3 then \code{plot.gbm} produces the plots. Otherwise, \code{plot.GBMFit}
#' returns only the grid of evaluation points and their average predictions
#' @param num_trees the number of trees used to generate the plot. Only the first
#' \code{num_trees} trees will be used
#' @param continuous_resolution The number of equally space points at which to
#' evaluate continuous predictors
#' @param grid_levels A list containing the points at which to evaluate each
#' predictor in \code{var_index} (in the same order as \code{var_index}). For
#' continuous predictors this is usually a regular sequence of values within
#' the range of the variable. For categorical predictors, the points are the
#' levels of the factor. When \code{length(var_index)} is one, the values can be
#' provided directly, outside a list.  This is NULL by default and generated
#' automatically from the data, using \code{continuous.resolution} for
#' continuous predictors. Forcing the values can be useful to evaluate two
#' models on the same exact range
#' @param return_grid if \code{TRUE} then \code{plot.GBMFit} produces no graphics
#' and only returns the grid of evaluation points and their average
#' predictions. This is useful for customizing the graphics for special
#' variable types or for dimensions greater than 3
#' @param type the type of prediction to plot on the vertical axis. See
#' \code{predict.GBMFit}
#' @param \dots other arguments passed to the plot function
#' @return Nothing unless \code{return_grid} is true then \code{plot.GBMFit}
#' produces no graphics and only returns the grid of evaluation points and
#' their average predictions.
#' @seealso \code{\link{gbm2}},
#' \code{\link[graphics]{plot}}
#' @references J.H. Friedman (2001). "Greedy Function Approximation: A Gradient
#' Boosting Machine," Annals of Statistics 29(4).
#' @keywords hplot
#' @export
plot.GBMFit <- function(gbm_fit_obj,
                        var_index=1,
                        n.trees=gbm_fit_obj$params$num_trees,
                        continuous_resolution=100,
                        grid_levels=NULL,
                        return_grid=FALSE,
                        type="link",
                        ...)
{
  # Check inputs
  check_if_gbm_fit(gbm_fit_obj)
  if (!is.element(type, c("link", "response"))){
    stop( "type must be either 'link' or 'response'")
  }
  
  if(all(is.character(var_index)))
  {
    i <- match(var_index, gbm_fit_obj$variables$var_names)
    if(any(is.na(i)))
    {
      stop("Plot variables not used in gbm model fit: ",i.var[is.na(i)])
    } else
    {
      var_index <- i
    }
  }
  
  if((min(var_index)<1) || (max(var_index) > length(gbm_fit_obj$variables$var.names)))
  {
    warning("var_index must be between 1 and ", length(gbm_fit_obj$variables$var.names))
  }
  if(num_trees > gbm_fit_obj$variables$n.trees)
  {
    warning(paste("num_trees exceeds the number of trees in the model, ", 
                  gbm_fit_obj$variables$num_trees,
                  ". Plotting using ", gbm_fit_obj$variables$n.trees," trees.",sep=""))
    num_trees <- gbm_fit_obj$variables$num_trees
  }
  
  if(length(var_index) > 3)
  {
    warning("gbm.int.plot creates up to 3-way interaction plots.\nplot.gbm will only return the plotting data structure.")
    return_grid = TRUE
  }
  
  # generate grid to evaluate gbm model
  if (is.null(grid_levels)) {
    grid_levels <- vector("list",length(var_index))
    for(i in seq_len(length(var_index)))
    {
      # continuous
      if(is.numeric(gbm_fit_obj$variables$var_levels[[var_index[i]]]))
      {
        grid_levels[[i]] <- seq(min(gbm_fit_obj$variables$var_levels[[var_index[i]]]),
                                max(gbm_fit_obj$variables$var_levels[[var_index[i]]]),
                                length=continuous_resolution)
      }
      # categorical or ordered
      else
      {
        grid_levels[[i]] <- as.numeric(factor(gbm_fit_obj$variables$var_levels[[var_index[i]]],
                                              levels=gbm_fit_obj$variables$var_levels[[var_index[i]]]))-1
      }
    }
  }
  else
  {
    # allow grid.levels to not be a list when there is only one predictor
    if (length(var_index) == 1 & !is.list(grid_levels)) {
      grid_levels <- list(grid_levels)
    }
    # check compatibility in length, at least
    if (length(grid_levels) != length(var_index)) {
      stop("Need grid_levels for all variables in var_index")
    }
    # convert levels for categorical predictors into numbers
    for(i in seq_len(length(var_index))) {
      if(!is.numeric(gbm_fit_obj$variables$var_levels[[var_index[i]]]))
      {
        grid_levels[[i]] <- as.numeric(grid_levels[[i]]) - 1
      }
    }
  }
  
  X <- expand.grid(grid_levels)
  names(X) <- paste("X", seq_len(length(var_index)),sep="")
  
  # Next if block for compatibility with objects created with 1.6
  if (is.null(gbm_fit_obj$num.classes)){
    gbm_fit_obj$num.classes <- 1
  }
  
  # evaluate at each data point
  y <- .Call("gbm_plot",
             X = data.matrix(X),
             i.var = as.integer(var_index-1),
             n.trees = as.integer(num_trees) ,
             initF = as.double(gbm_fit_obj$initF),
             trees = gbm_fit_obj$trees,
             c.splits = gbm_fit_obj$c.splits,
             var.type = as.integer(gbm_fit_obj$variables$var_type),
             PACKAGE = "gbm")
  
  if(is.element(gbm_fit_obj$distribution$name, c("bernoulli", "pairwise")) && type=="response") {
    X$y <- 1/(1+exp(-y))
  }
  else if ((gbm_fit_obj$distribution$name=="poisson") && (type=="response")){
    X$y <- exp(y)
  }
  else if ((gbm_fit_obj$distribution$name=="gamma") && (type=="response")){
    X$y <- exp(y)
  }
  else if ((gbm_fit_obj$distribution$name=="tweedie") && (type=="response")){
    X$y <- exp(y)
  }
  else if (type=="response"){
    warning("type 'response' only implemented for 'bernoulli', 'poisson', 'gamma', 'tweedie', and 'pairwise'. Ignoring" )
  }
  else { X$y <- y }
  
  # transform categorical variables back to factors
  f.factor <- rep(FALSE,length(var_index))
  for(i in seq_len(length(var_index)))
  {
    if(!is.numeric(gbm_fit_obj$variables$var_levels[[var_index[i]]]))
    {
      X[,i] <- factor(gbm_fit_obj$variables$var_levels[[var_index[i]]][X[,i]+1],
                      levels=gbm_fit_obj$variables$var_levels[[var_index[i]]])
      f.factor[i] <- TRUE
    }
  }
  
  if(return_grid)
  {
    names(X)[seq_len(length(var_index))] <- gbm_fit_obj$variables$var_names[var_index]
    return(X)
  }
  
  # create the plots
  if(length(var_index)==1)
  {
    if(!f.factor)
    {
      j <- order(X$X1)
      
      if (is.element(gbm_fit_obj$distribution$name, c("bernoulli", "pairwise"))) {
        if ( type == "response" ){
          ylabel <- "Predicted probability"
        }
        else {
          ylabel <- paste("f(", gbm_fit_obj$variables$var_names[var_index],")",sep="")
        }
        plot( X$X1, X$y , type = "l", xlab = gbm_fit_obj$variables$var_names[var_index], ylab=ylabel )
      }
      else if ( gbm_fit_obj$distribution$name == "poisson" ){
        if (type == "response" ){
          ylabel <- "Predicted count"
        }
        else{
          ylabel <- paste("f(",gbm_fit_obj$variables$var_names[var_index],")",sep="")
        }
        plot( X$X1, X$y , type = "l", xlab = gbm_fit_obj$variables$var_names[var_index], ylab=ylabel )
      }
      else {
        plot(X$X1,X$y,
             type="l",
             xlab=gbm_fit_obj$variables$var_names[var_index],
             ylab=paste("f(", gbm_fit_obj$variables$var_names[var_index],")",sep=""),...)
      }
    }
    else
    {
      if (is.element(gbm_fit_obj$distribution$name, c("bernoulli", "pairwise")) && type == "response" ){
        ylabel <- "Predicted probability"
        plot( X$X1, X$y, type = "l", xlab=gbm_fit_obj$variables$var_names[var_index], ylab=ylabel )
      }
      else if ( gbm_fit_obj$distribution$name == "poisson" & type == "response" ){
        ylabel <- "Predicted count"
        plot( X$X1, X$y, type = "l", xlab=gbm_fit_obj$variables$var_names[var_index], ylab=ylabel )
      }
      else {
        plot(X$X1,X$y,
             type="l",
             xlab=gbm_fit_obj$variables$var_names[var_index],
             ylab=paste("f(", gbm_fit_obj$variables$var_names[var_index],")",sep=""),...)
      }
    }
  }
  else if(length(var_index)==2)
  {
    if(!f.factor[1] && !f.factor[2])
    {
      print(levelplot(y~X1*X2,data=X,
                      xlab=gbm_fit_obj$variables$var_names[var_index[1]],
                      ylab=gbm_fit_obj$variables$var.names[var_index[2]],...))
      
    }
    else if(f.factor[1] && !f.factor[2])
    {
      print(xyplot(y~X2|X1,data=X,
                   xlab=gbm_fit_obj$variables$var_names[var_index[2]],
                   ylab=paste("f(", gbm_fit_obj$variables$var_names[var_index[1]],",", gbm_fit_obj$variables$var_names[var_index[2]],")",sep=""),
                   type="l",
                   panel = panel.xyplot,
                   ...))
    }
    else if(!f.factor[1] && f.factor[2])
    {
      print(xyplot(y~X1|X2,data=X,
                   xlab=gbm_fit_obj$variables$var_names[var_index[1]],
                   ylab=paste("f(",gbm_fit_obj$variables$var_names[var_index[1]],",",gbm_fit_obj$variables$var_names[var_index[2]],")",sep=""),
                   type="l",
                   panel = panel.xyplot,
                   ...))
    }
    else
    {
      print(stripplot(X1~y|X2,data=X,
                      xlab=gbm_fit_obj$variables$var_names[var[2]],
                      ylab=paste("f(",gbm_fit_obj$variables$var_names[var_index[1]],",",gbm_fit_obj$variables$var_names[var_index[2]],")",sep=""),
                      ...))
      
    }
  }
  else if(length(var_index)==3)
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
                      xlab=gbm_fit_obj$variables$var_names[var_index[i[1]]],
                      ylab=gbm_fit_obj$variables$var_names[var_index[i[2]]],...))
      
    }
    # 1 factor, 2 continuous
    else if(sum(f.factor)==1)
    {
      print(levelplot(y~X1*X2|X3,data=X.new,
                      xlab=gbm_fit_obj$variables$var_names[var_index[i[1]]],
                      ylab=gbm_fit_obj$variables$var_names[var_index[i[2]]],...))
      
    }
    # 2 factors, 1 continuous
    else if(sum(f.factor)==2)
    {
      print(xyplot(y~X1|X2*X3,data=X.new,
                   type="l",
                   xlab=gbm_fit_obj$variables$var_names[var_index[i[1]]],
                   ylab=paste("f(",paste(gbm_fit_obj$variables$var_names[var_index[1:3]],collapse=","),")",sep=""),
                   panel = panel.xyplot,
                   ...))
    }
    # 3 factors, 0 continuous
    else if(sum(f.factor)==3)
    {
      print(stripplot(X1~y|X2*X3,data=X.new,
                      xlab=gbm_fit_obj$variables$var_names[var_index[i[1]]],
                      ylab=paste("f(",paste(gbm_fit_obj$variables$var_names[var_index[1:3]],collapse=","),")",sep=""),
                      ...))
      
    }
  }
}
