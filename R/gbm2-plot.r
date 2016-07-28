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
#' @export plot.GBMFit
plot.GBMFit <- function(gbm_fit_obj,
                        var_index=1,
                        num_trees=gbm_fit_obj$params$num_trees,
                        continuous_resolution=100,
                        grid_levels=NULL,
                        return_grid=FALSE,
                        type="link",
                        ...)
{
  # Check inputs
  if (!is.element(type, c("link", "response"))){
    stop( "type must be either 'link' or 'response'")
  }
  
  if(all(is.character(var_index))) {
    i <- match(var_index, gbm_fit_obj$variables$var_names)
    if(any(is.na(i))) {
      stop("Plot variables not used in gbm model fit: ",i.var[is.na(i)])
    } else {
      var_index <- i
    }
  }
  
  if((min(var_index)<1) || (max(var_index) > length(gbm_fit_obj$variables$var_names))) {
    warning("var_index must be between 1 and ", length(gbm_fit_obj$variables$var_names))
  }
  if(num_trees > gbm_fit_obj$variables$num_trees) {
    warning(paste("num_trees exceeds the number of trees in the model, ", 
                  gbm_fit_obj$variables$num_trees,
                  ". Plotting using ", gbm_fit_obj$variables$num_trees," trees.",sep=""))
    num_trees <- gbm_fit_obj$variables$num_trees
  }
  
  if(length(var_index) > 3) {
    warning("gbm.int.plot creates up to 3-way interaction plots.\nplot.gbm will only return the plotting data structure.")
    return_grid = TRUE
  }
  
  # generate grid to evaluate gbm model
  if (is.null(grid_levels)) {
   grid_levels <- get_default_grid_levels(gbm_fit_obj, var_index, continuous_resolution)
  } else {
    grid_levels <- generate_grid_levels(grid_levels, gbm_fit_obj, var_index)
  }
  
  # Expand the grid of variables
  X <- expand.grid(grid_levels)
  names(X) <- paste("X", seq_len(length(var_index)),sep="")

  # evaluate at each data point
  y <- .Call("gbm_plot",
             X = data.matrix(X),
             i.var = as.integer(var_index-1),
             num_trees = as.integer(num_trees) ,
             initF = as.double(gbm_fit_obj$initF),
             trees = gbm_fit_obj$trees,
             c.splits = gbm_fit_obj$c.splits,
             var.type = as.integer(gbm_fit_obj$variables$var_type),
             PACKAGE = "gbm")
  
  if(type=="response") {
   X$y <- response(y, gbm_fit_obj$distribution)
  } else { 
    X$y <- y 
  }
  
  # transform categorical variables back to factors
  f.factor <- rep(FALSE, length(var_index))
  for(i in seq_len(length(var_index))) {
    if(!is.numeric(gbm_fit_obj$variables$var_levels[[var_index[i]]])) {
      X[,i] <- factor(gbm_fit_obj$variables$var_levels[[var_index[i]]][X[,i]+1],
                      levels=gbm_fit_obj$variables$var_levels[[var_index[i]]])
      f.factor[i] <- TRUE
    }
  }
  
  # Stop if grid is returned
  if(return_grid) {
    names(X)[seq_len(length(var_index))] <- gbm_fit_obj$variables$var_names[var_index]
    return(X)
  }
  
  # create the plots
  if(length(var_index)==1) {
      ylabel <- paste("f(", gbm_fit_obj$variables$var_names[var_index],")",sep="")
      if(type=="response") ylabel <- get_ylabel_one_var(gbm_fit_obj$distribution)
      plot(X$X1,X$y,
          type="l",
          xlab=gbm_fit_obj$variables$var_names[var_index],
          ylab=ylabel, ...)
 
  } else if(length(var_index)==2) {
    select_two_var_plot(f.factor, X, gbm_fit_obj, var_index, ...)
  } else if(length(var_index)==3) {
    select_three_var_plot(f.factor, X, gbm_fit_obj, var_index, ...)
  }
}


##### Helper Functions #####
#' @export
get_default_grid_levels <- function(gbm_fit_obj, var_index, continuous_resolution) {
  # Vector of lists is output
  grid_levels <- vector("list",length(var_index))
  
  # Loop over variables and get levels
  for(i in seq_len(length(var_index))) {
    # continuous
    if(is.numeric(gbm_fit_obj$variables$var_levels[[var_index[i]]])) {
      grid_levels[[i]] <- seq(min(gbm_fit_obj$variables$var_levels[[var_index[i]]]),
                              max(gbm_fit_obj$variables$var_levels[[var_index[i]]]),
                              length=continuous_resolution)
    }
    # categorical or ordered - convert to appropriate numerics
    else {
      grid_levels[[i]] <- as.numeric(factor(gbm_fit_obj$variables$var_levels[[var_index[i]]],
                                            levels=gbm_fit_obj$variables$var_levels[[var_index[i]]]))-1
    }
  }
  
  return(grid_levels)
}

#' @export
generate_grid_levels <- function(grid_levels, gbm_fit_obj, var_index) {
  
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
    if(!is.numeric(gbm_fit_obj$variables$var_levels[[var_index[i]]])) {
      grid_levels[[i]] <- as.numeric(grid_levels[[i]]) - 1
    }
  }
  
  return(grid_levels)
}

#' @export response
response <- function(resp, dist_obj) {
  UseMethod("response", dist_obj)
}

response.default <- function(resp, dist_obj) {
  warning("type 'response' only implemented for 'Bernoulli', 'Poisson', 'Gamma', 'Tweedie', and 'Pairwise'. Ignoring" )
  return(NULL)
}

response.BernoulliGBMDist <- function(resp, dist_obj) {
  return(1/(1+exp(-resp)))
}

response.GammaGBMDist <- function(resp, dist_obj) {
  return(exp(resp))
}

response.PairwiseGBMDist <- function(resp, dist_obj) {
  return(1/(1+exp(-resp)))
}

response.PoissonGBMDist <- function(resp, dist_obj) {
  return(exp(resp))
}

response.TweedieGBMDist <- function(resp, dist_obj) {
  return(exp(resp))
}