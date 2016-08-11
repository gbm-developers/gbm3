# Create Strata 
# 
# Function used to create to appropriate strata vectors -
# currently only implemented for the CoxPH distribution.
# 
# @usage create_strata(gbm_data_obj, train_params, distribution_obj,
#  order_strata_by_id=TRUE)
# 
# @param gbm_data_obj gbm data object
# 
# @param train_params gbm training parameters object
# 
# @param distribution_obj a gbm distribution object - the strata in this object will be
# updated if they exist - this is only the case for the CoxPH model at this moment.
# 
# @param order_strata_by_id indicates whether or not to order the entire strata
# vector by observation id. This defaults to \code{TRUE} but should be set to \code{FALSE}
# for cross-validation strata creation. 
# 
# @author James Hickey
# 
# @return updated distribution_obj
# 

create_strata <- function(gbm_data_obj, train_params, distribution_obj, order_strata_by_id=TRUE) {
  check_if_gbm_dist(distribution_obj)
  check_if_gbm_data(gbm_data_obj)
  
  # Put in defaults
  if(is.null(distribution_obj$sorted)) distribution_obj$sorted <- NA
  if(is.null(distribution_obj$strata)) distribution_obj$strata <- NA
  
  if(distribution_obj$name == "CoxPH") {
    num_train_rows <- train_params$num_train_rows
    num_test_rows <- nrow(gbm_data_obj$x) - num_train_rows
    
    # Determine test indices
    if(num_test_rows == 0) {
      test_indices <- 0
    } else {
      test_indices <- (num_train_rows+1):nrow(gbm_data_obj$x)
    }
      
    # Set up strata 
    if(!is.na(distribution_obj$original_strata_id[1])) {
      # Sort strata according to patient ID
      distribution_obj$strata <- distribution_obj$original_strata_id
      if(order_strata_by_id) {
        distribution_obj$strata <- distribution_obj$strata[order(train_params$id)]
      }
      
      # Order strata and split into train/test
      strataVecTrain <- distribution_obj$strata[seq_len(num_train_rows)]
      strataVecTest <- distribution_obj$strata[test_indices]
      
      # Cum sum the number in each stratum and pad with NAs
      # between train and test strata
      strataVecTrain <- as.vector(cumsum(table(strataVecTrain)))
      strataVecTest <- as.vector(cumsum(table(strataVecTest)))
      
      strataVecTrain <- c(strataVecTrain, rep(NA, num_train_rows-length(strataVecTrain)))
      
      # If no test set make empty
      if(num_test_rows == 0) {
        strataVecTest <- c() 
      } else {
        strataVecTest <- c(strataVecTest, rep(NA, max(num_test_rows-length(strataVecTest), 0)))  
      }
    
      # Recreate Strata Vec to Pass In
      nstrat <- c(strataVecTrain, strataVecTest)
      
    }
    else
    {
      # Put all the train and test data in a single stratum
      distribution_obj$strata <- rep(1, nrow(gbm_data_obj$x))
      trainStrat <- c(num_train_rows, rep(NA, num_train_rows-1))
      if(num_test_rows == 0) {
        testStrat <- c()
      } else {
        testStrat <- c(num_test_rows, rep(NA, max(num_test_rows-1, 0)))
      }
      nstrat <- c(trainStrat, testStrat)
    }
    
    # Sort response according to strata
    # i_order sets order of outputs
    if (ncol(gbm_data_obj$y)==2) {
      sorted <- c(order(distribution_obj$strata[seq_len(num_train_rows)], -gbm_data_obj$y[seq_len(num_train_rows), 1]),
                  order(distribution_obj$strata[test_indices],
                        -gbm_data_obj$y[test_indices, 1])) 
      i_order <- c(order(distribution_obj$strata[seq_len(num_train_rows)], -gbm_data_obj$y[1:num_train_rows, 1]),
                   order(distribution_obj$strata[test_indices],
                         -gbm_data_obj$y[test_indices, 1]) + num_train_rows)
    } else if (ncol(gbm_data_obj$y)==3) {
      sorted <- cbind(c(order(distribution_obj$strata[seq_len(num_train_rows)], -gbm_data_obj$y[seq_len(num_train_rows), 1]),
                        order(distribution_obj$strata[test_indices], -gbm_data_obj$y[test_indices, 1])),
                      c(order(distribution_obj$strata[seq_len(num_train_rows)], -gbm_data_obj$y[seq_len(num_train_rows), 2]),
                        order(distribution_obj$strata[test_indices], -gbm_data_obj$y[test_indices, 2])))
      i_order <- c(order(distribution_obj$strata[seq_len(num_train_rows)], -gbm_data_obj$y[seq_len(num_train_rows), 1]),
                   order(distribution_obj$strata[test_indices], -gbm_data_obj$y[test_indices, 1]) + num_train_rows)
    }
    else
    {
      stop("Survival object must be either right or counting type.")
    }
    
    
    # Add in sorted column and strata
    StrataVec <-  nstrat
    sortedVec <- sorted-1L
    
    distribution_obj$time_order <- i_order
    distribution_obj$sorted <- as.matrix(as.data.frame(sortedVec))
    distribution_obj$strata <- as.double(StrataVec)
    
  }
  
  return(distribution_obj)
}

