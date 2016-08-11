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