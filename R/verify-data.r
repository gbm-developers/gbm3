# Verify data
# 
# Functions used to determine whether the data input is valid - checks that
# the predictor and response vectors are of the right form
# 
# @usage verify-data(x, y)
# 
# @param x the matrix or data-frame containing the predictive variables
# 
# @param y the response vector
#
# @author James Hickey
#
# @return NA - Throws an error if data is of incorrect type
# 

verify_data <- function(x, y) {
  # Check predictor variables are either data-frame or matrix
  if(!is.matrix(x) && !is.data.frame(x)) {
    stop("Predictor variables, x, must be contained in a dataframe or matrix")
  }
  
  # Check same for y
  if(!is.matrix(y) && !is.data.frame(y) && !is.atomic(y)) {
    stop("Response, y, must be contained in a dataframe, matrix or vector")
  }
  
  # Check same length
  check_sanity(x, y)
  
  # Check variables are either numeric, ordered or factor
  check_var_type(x, y)
  
  # Check missing values are of correct form
  checkMissing(x, y)
  
}