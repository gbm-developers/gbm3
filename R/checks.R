checkMissing <- function(x, y){
   nms <- getVarNames(x)
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

checkID <- function(id){
   # Check for disallowed interaction.depth
   if(id < 1) {
      stop("interaction.depth must be at least 1.")
   }
   else if(id > 49) {
      stop("interaction.depth must be less than 50. You should also ask yourself why you want such large interaction terms. A value between 1 and 5 should be sufficient for most applications.")
   }
   invisible(id)
}

checkWeights <- function(w, n){
   # Logical checks on weights
   if(length(w)==0) { w <- rep(1, n) }
   else if(any(w < 0)) stop("negative weights not allowed")
   w
}

checkOffset <- function(o, y, dist){
   # Check offset
  if(is.null(o))
      o <- rep(0,length(y))
   else if((length(o) != length(y)) & dist != "coxph")
      stop("The length of offset does not equal the length of y.")
   else if ((length(o) != (length(y)/2)) & dist == "coxph")
     stop("The length of offset does not equal the length of y.")
   else if(!is.numeric(o))
     stop("offset must be numeric")
   else if(sum(is.na(o))>0)
     stop("offset can not contain NA's")

   o
}

getVarNames <- function(x){
   if(is.matrix(x)) { var.names <- colnames(x) }
   else if(is.data.frame(x)) { var.names <- names(x) }
   else { var.names <- paste("X",1:ncol(x),sep="") }
   var.names
}



checkSanity <- function(x, y){
  
  nms <- getVarNames(x)
  
  # x and y are not the same length
  if(nrow(x) != ifelse("Surv" %in% class(y), nrow(y), length(y))) {
    stop("The number of rows in x does not equal the length of y.")
  }
  
}

checkVarType <- function(x, y){
  
  nms <- getVarNames(x)
  
  # Excessive Factors
  Factors <- sapply(x, function(X){is.factor(X)})
  nLevels <- sapply(x, function(X){length(levels(X))})
  
  excessLevels <- nLevels > 1024
  excessLevelsIndex <- paste(which(excessLevels), collapse = ', ')
  excessLevelsVars <- paste(nms[which(excessLevels)], collapse = ', ')
  
  if(any(excessLevels)) {
    stop("gbm does not currently handle categorical variables with more than 1024 levels. Variable ", excessLevelsIndex,": ", excessLevelsVars," has ", nLevels[which(excessLevels)]," levels.")
    
  }
  
  # Not an acceptable class
  inacceptClass <- sapply(x, function(X){! (is.ordered(X) | is.factor(X) | is.numeric(X)) })
  inacceptClassIndex <- paste(which(inacceptClass), collapse = ', ')
  inacceptClassVars <- paste(nms[which(inacceptClass)], collapse = ', ')
  
  if(any(inacceptClass)){
    stop("variable ", inacceptClassIndex,": ", inacceptClassVars, " is not of type numeric, ordered, or factor.")
  }
  
}

