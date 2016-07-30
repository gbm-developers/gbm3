# Series of internal functions used 
# to check inputs and convert parameters

##### Check gbm objects #####
check_if_gbm_dist <- function(distribution_obj) {
  # Check if GBM dist object
  if(!any(class(distribution_obj) %in% paste0(available_distributions(), "GBMDist"))) {
    stop("Function requires a specific GBMDist object.")
  }
}

check_if_gbm_data <- function(data_obj) {
  if(!any(class(data_obj) %in% "GBMData")) {
    stop("Function requires a GBMData object.")
  }
}

check_if_gbm_fit <- function(fit_obj) {
  if(!any(class(fit_obj) %in% "GBMFit")) {
    stop("Function requires a GBMFit object.")
  }
}

check_if_gbm_train_params <- function(params_obj) {
  if(!any(class(params_obj) %in% "GBMTrainParams")) {
    stop("Function requires a GBMTrainParams object.")
  }
}

check_if_gbm_var_container <- function(var_obj) {
  if(!any(class(var_obj) %in% "GBMVarCont")) {
    stop("Function requires a GBMVarCont object.")
  }
}

#### Check function inputs ####
check_cv_parameters <- function(cv_folds, cv_class_stratify, fold_id, train_params) {
  check_if_natural_number(cv_folds, "cv_folds")
  if(!is.logical(cv_class_stratify)) stop("cv_class_stratify must be a logical")
  check_if_gbm_train_params(train_params)

  # Check fold_id does not split observation data up
  if(!is.null(fold_id)) {
    for(id in train_params$id) {
      if(length(unique(fold_id[train_params$id == id])) > 1) {
        stop("Observations are split across multiple folds")
      }
    }
  }
}

check_if_natural_number <- function(value, name) {
  # value - the value of the parameter to check 
  if(is.null(value) || is.infinite(value)  || is.logical(value)
     || !(abs(value - round(value)) < .Machine$double.eps^0.5) ||
     (value <= 0) || (length(value) > 1)) {
    stop("The parameter  ", name,  " must be a positive whole number")
  }
}

checkMissing <- function(x, y){
   nms <- get_var_names(x)
   #### Check for NaNs in x and NAs in response
   j <- apply(x, 2, function(z) any(is.nan(z)))
   if(any(j)) {
      stop("Use NA for missing values. NaN found in predictor variables:",
           paste(nms[j],collapse=","))
   }
   
   if(any(is.na(y))) stop("Missing values are not allowed in the response")
   
   AllMiss <- apply(x, 2, function(X){all(is.na(X))})
   AllMissVarIndex <- paste(which(AllMiss), collapse = ', ')
   AllMissVar <- paste(nms[which(AllMiss)], collapse = ', ')
   
   if(any(AllMiss)) {
      stop("variable(s) ", AllMissVarIndex, ": ", AllMissVar, " contain only missing values.")
   }
   
   invisible(NULL)
 }

check_interaction_depth <- function(id){
   # Check for disallowed interaction.depth
   if(id < 1) {
      stop("interaction_depth must be at least 1.")
   }
   else if(id > 49) {
      stop("interaction_depth must be less than 50. You should also ask yourself why you want such large interaction terms. A value between 1 and 5 should be sufficient for most applications.")
   }
   invisible(id)
}

check_weights <- function(w, n){
   # Logical checks on weights
   if(length(w)==0) { w <- rep(1, n) }
   else if(any(w != as.double(w))) {stop("weights must be doubles")}
   else if(any(w < 0)) {stop("negative weights not allowed")}
   w
}

check_offset <- function(o, y, dist){
   # Check offset
  if(is.null(o))
      o <- rep(0,length(y))
   else if((length(o) != length(y)) && (dist$name != "CoxPH"))
      stop("The length of offset does not equal the length of y.")
   else if(!is.numeric(o))
     stop("offset must be numeric")
   else if(sum(is.na(o))>0)
     stop("offset can not contain NA's")

   o
}

get_var_names <- function(x){
  if(is.matrix(x)) { var.names <- colnames(x) }
  else if(is.data.frame(x)) { var.names <- names(x) }
  else { var.names <- paste("X", 1:ncol(x),sep="") }
  var.names
}

check_sanity <- function(x, y){
  # x and y are not the same length
  if(nrow(x) != ifelse(!is.null(dim(y)), nrow(y), length(y))) {
    stop("The number of rows in x does not equal the length of y.")
  }
}

check_var_type <- function(x, y){
  
  nms <- get_var_names(x)
  
  # Excessive Factors
  Factors <- vapply(x, is.factor, TRUE)
  nLevels <- vapply(x, nlevels, 0L)
  
  excessLevels <- nLevels > 1024
  excessLevelsIndex <- paste(which(excessLevels), collapse = ', ')
  excessLevelsVars <- paste(nms[which(excessLevels)], collapse = ', ')
  
  if(any(excessLevels)) {
    stop("gbm does not currently handle categorical variables with more than 1024 levels. Variable ", excessLevelsIndex,": ", excessLevelsVars," has ", nLevels[which(excessLevels)]," levels.")
    
  }
  
  # Not an acceptable class
  inacceptClass <- vapply(x, function(X){! (is.ordered(X) | is.factor(X) | is.numeric(X)) }, TRUE)
  inacceptClassIndex <- paste(which(inacceptClass), collapse = ', ')
  inacceptClassVars <- paste(nms[which(inacceptClass)], collapse = ', ')
  
  if(any(inacceptClass)){
    stop("variable ", inacceptClassIndex,": ", inacceptClassVars, " is not of type - numeric, ordered or factor.")
  }
  
}

#### Data conversion functions ####
convertY <- function(y){

  FactorsY <- is.factor(y)
  nLevelsY <- nlevels(y)
  
  if(FactorsY & nLevelsY == 2){
    Y = as.numeric(y == levels(y)[2])
  } else {
    Y = y
  }
  
  return(Y)
  
}

## miscellaneous small internal helper functions.

## these are not exported and not formally documented in the manual

## Warn if a variable does not vary
##
## @param x a numeric variable to check
## @param ind an index to use in a potential warning
## @param name a name to include in a potential warning
## @return NULL, invisibly

warnNoVariation <- function(x, ind, name) {
  ## suppress warnings in here
  ## because min and max warn if the variable is completely NA
  suppressWarnings(variation <- range(x, na.rm=TRUE))
  
  ## I really mean ">=" here, which catches the all NA case
  ## and the standard case
  if (variation[[1]] >= variation[[2]]) {
    warning("variable ", ind, ": ", name, " has no variation.")
  }
  
  invisible(NULL)
}

convert_strata <- function(strata) {
  # If factor then convert to integer
  if(is.factor(strata)) {
    strata <- as.integer(strata)
  }
  
  # If it isn't default then check
  if(!is.na(strata[1])) {
    if(!is.vector(strata) || any(is.character(strata)) || any(is.infinite(strata)) || any(is.nan(strata)) ||
       !(all(strata == as.factor(strata)) || all(strata == as.integer(strata)))) {
      stop("strata must be an atomic vector of factors or integers")
    }
  }
  
  return(strata)
}

guess_distribution <- function(response) {
  # This function guesses the distribution if one is not provided
  if(length(unique(response)) == 2) {
    name <- "Bernoulli"
  } else if (class(response) == "Surv") {
    name <- "CoxPH"
  } else {
    name <- "Gaussian"
  }
  message(paste("Distribution not specified, assuming", name, "...\n"))
  return(list(name=name))
}

to_old_gbm <- function(gbm_fit_obj) {
  # Check input
  check_if_gbm_fit(gbm_fit_obj)
  
  # Convert gbm_fit_obj to old API
  gbm_fit_old <- list()
  gbm_fit_old$initF <- gbm_fit_obj$initF
  gbm_fit_old$train.error <- gbm_fit_obj$train.error
  gbm_fit_old$valid.error <- gbm_fit_obj$valid.error
  gbm_fit_old$trees <- gbm_fit_obj$trees
  gbm_fit_old$c.splits <- gbm_fit_obj$c.splits
  gbm_fit_old$oobag.improve <- gbm_fit_obj$oobag.improve
  gbm_fit_old$fit <- gbm_fit_obj$fit
  
  gbm_fit_old$bag.fraction <- gbm_fit_obj$params$bag_fraction
  gbm_fit_old$distribution <- tolower(gbm_fit_obj$distribution$name)
  gbm_fit_old$interaction.depth <- gbm_fit_obj$params$interaction_depth
  gbm_fit_old$n.minobsinnode <- gbm_fit_obj$params$min_num_obs_in_node
  gbm_fit_old$n.trees <- length(gbm_fit_obj$trees)
  gbm_fit_old$nTrain <- gbm_fit_obj$params$num_train_rows
  gbm_fit_old$nTrainPats <- gbm_fit_obj$params$num_train
  gbm_fit_old$patient.id <- gbm_fit_obj$params$id
  gbm_fit_old$mFeatures <- gbm_fit_obj$params$num_features
  gbm_fit_old$train.fraction <- gbm_fit_obj$params$train_fraction
  gbm_fit_old$response.name <- "Response" # Default value
  gbm_fit_old$shrinkage <- gbm_fit_obj$params$shrinkage
  gbm_fit_old$var.levels <- gbm_fit_obj$variables$var_levels
  gbm_fit_old$var.monotone <- gbm_fit_obj$variables$var_monotone
  gbm_fit_old$var.names <- gbm_fit_obj$variables$var_names
  gbm_fit_old$var.type <- gbm_fit_obj$variables$var_type
  gbm_fit_old$verbose <- FALSE # Default value
  gbm_fit_old$strata <- gbm_fit_obj$distribution$strata
  gbm_fit_old$sorted <- gbm_fit_obj$distribution$sorted
  gbm_fit_old$prior.node.coeff.var <- gbm_fit_obj$distribution$prior_node_coeff_var
  
  if(!is.null(gbm_fit_obj$gbm_data_obj)) {
    # put the observations back in - these are ordered according to id and group
    data <- gbm_fit_obj$gbm_data_obj
    gbm_fit_old$data <- list(y=data$y, x=data$x,x.order=data$x_order, offset=data$offset, 
                             Misc=unlist(get_misc(gbm_fit_obj$distribution)), w=data$weights)
  } else {
    gbm_fit_old$data <- NULL
  }
  
  gbm_fit_old$cv.folds <- gbm_fit_obj$cv_folds
  gbm_fit_old$cv.error <- gbm_fit_obj$cv_error
  gbm_fit_old$cv.fitted <- gbm_fit_obj$cv.fitted
  gbm_fit_old$Terms <- gbm_fit_obj$Terms
  gbm_fit_old$call <- gbm_fit_obj$call
  gbm_fit_old$m <- gbm_fit_obj$m
  class(gbm_fit_old) <- "gbm"
  
  message("Converted to old gbm object - this will not work with new packages functions")
  
  return(gbm_fit_old)
  
}