#' Verify data
#' 
#' Functions used to determine whether the data input is valid - checks that
#' the predictor and response vectors are of the right form
#' 
#' @usage verify-data(x, y)
#' 
#' @param x the matrix or data-frame containing the predictive variables
#' 
#' @param y the response vector
#' 
#' @return NA - Throws an error if data is of incorrect type
#' 

verify_data <- function(x, y) {
  # Check predictor variables are either data-frame or matrix
  if(!is.matrix(x) || !is.data.frame(x)) {
    stop("Predictor variables, x, must be contained in a dataframe or matrix")
  }
  
  # Check same for y
  if(!is.matrix(y) || !is.data.frame(y) || !is.atomic(y)) {
    stop("Response, x, must be contained in a dataframe, matrix or vector")
  }
  
  # Check same length
  checkSanity(x, y)
  
  # Check variables are either numeric, ordered or factor
  checkVarType(x, y)
  
  # Check missing values are of correct form
  checkMissing(x, y)
  
}