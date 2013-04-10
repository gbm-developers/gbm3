checkMissing <- function(x, y){
   nms <- getVarNames(x)
   #### Check for NaNs in x and NAs in response
   j <- apply(x, 2, function(z) any(is.nan(z)))
   if(any(j)) {
      stop("Use NA for missing values. NaN found in predictor variables:",
           paste(nms[j],collapse=","))
   }
   if(any(is.na(y))) stop("Missing values are not allowed in the response")
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

checkOffset <- function(o, y){
   # Check offset
   if(is.null(o) | all(o==0)) { o <- NA  }
   else if(length(o) != length(y))   {
      stop("The length of offset does not equal the length of y.")
   }
   o
}

getVarNames <- function(x){
   if(is.matrix(x)) { var.names <- colnames(x) }
   else if(is.data.frame(x)) { var.names <- names(x) }
   else { var.names <- paste("X",1:ncol(x),sep="") }
   var.names
 }
