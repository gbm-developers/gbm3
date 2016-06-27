#' Variable Container
#' 
#' Contains properties associated with predictor variables.
#' 
#' 
#' @usage var_container(gbm_data_obj,  var_monotone=NULL, var_names=NULL)
#' 
#' @param gbm_data_obj an initialized GBMData object for the fit
#'
#' @param var_monotone  an optional vector, the same length as the
#' number of predictors, indicating which variables have a monotone
#' increasing (+1), decreasing (-1), or arbitrary (0) relationship
#' with the outcome.
#' 
#' @param var_names a vector of strings of the same length as the 
#' 1st dimension of the response.
#' 
#' @return a GBMVarCont object
#' 
#' @export var_container
#' 

var_container <- function(gbm_data_obj, var_monotone=NULL, var_names=NULL) {
  check_if_gbm_data(gbm_data_obj)
  
  # Check var_monotone
  cCols <- ncol(gbm_data_obj$x)
  if(is.null(var_monotone)) var_monotone <- rep(0, cCols)
  else if(length(var_monotone)!=cCols)
  {
    stop("Length of var_monotone != number of predictors")
  }
  else if(!all(is.element(var_monotone,-1:1)))
  {
    stop("var_monotone must be -1, 0, or 1")
  }
  
  # Check names and Get
  if(is.null(var_names)) var_names <- getVarNames(gbm_data_obj$x)
  
  if(!is.null(var_names) && (!is.atomic(var_names) || any(var_names != as.character(var_names))
     || is.null(var_names)) ){
    stop("Names of data must be a vector of strings.")
  }
  
  if(!is.null(var_names) && length(var_names)!=cCols) stop("Length of var_names != number of predictors")
  
  # setup variable types
  var_type <- rep(0, cCols)
  var_levels <- vector("list", cCols)
  
  for(i in seq_len(var_type)) {
    if(is.ordered(gbm_data_obj$x[,i])) {
      
      var_levels[[i]] <- levels(factor(gbm_data_obj$x[,i]))
      var_type[i] <- 0
      
    } else if(is.factor(gbm_data_obj$x[,i])) {
      var_levels[[i]] <- levels(factor(gbm_data_obj$x[,i]))
      var_type[i] <- max(gbm_data_obj$x[,i],na.rm=TRUE)+1
      
    } else if(is.numeric(gbm_data_obj$x[,i])) {
      var_levels[[i]] <- quantile(gbm_data_obj$x[,i],prob=(0:10)/10,na.rm=TRUE)
    }
  }
  
  return(structure(list(var_monotone=var_monotone, var_names=var_names,
                        var_levels=var_levels, var_type=var_type),
                   class = "GBMVarCont"))
}

