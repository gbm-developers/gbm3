#' Create Strata 
#' 
#' Function used to create to appropriate strata vectors -
#' currently only implemented for the CoxPH distribution.
#' 
#' @usage create_strata(gbm_data_obj, train_params distribution_obj)
#' 
#' @param gbm_data_obj gbm data object
#' 
#' @param train_params gbm training parameters object
#' 
#' @param distribution_obj a gbm distribution object - the strata in this object will be
#' updated if they exist - this is only the case for the CoxPH model at this moment.
#' 
#' @return updated distribution_obj
#' 
#' @export

create_strata <- function(gbm_data_obj, train_params, distribution_obj) {
  check_if_gbm_dist(distribution_obj)
  check_if_gbm_data(gbm_data_obj)
  
  # Put in defaults
  distribution_obj$sorted <- NA
  distribution_obj$strata <- NA
  
  if(distribution_obj$name == "CoxPH") {
    num_train_rows <- sum(train_params$num_rows_per_obs[seq_len(train_params$num_train)])
    num_test_rows <- nrow(gbm_data_obj$x) - num_train_rows
      
    # Set up strata 
    if(!is.null(distribution_obj$strata))
    {
      # Sort strata according to patient ID
      distribution_obj$strata <- distribution_obj$strata[order(train_params$id)]
      
      # Order strata and split into train/test
      strataVecTrain <- distribution_obj$strata[1:num_train_rows]
      strataVecTest <- distribution_obj$strata[(num_train_rows+1): nrow(gbm_data_obj$x)]
      
      # Cum sum the number in each stratum and pad with NAs
      # between train and test strata
      strataVecTrain <- as.vector(cumsum(table(strataVecTrain)))
      strataVecTest <- as.vector(cumsum(table(strataVecTest)))
      
      strataVecTrain <- c(strataVecTrain, rep(NA, num_train_rows-length(strataVecTrain)))
      strataVecTest <- c(strataVecTest, rep(NA, num_test_rows-length(strataVecTest)))
      
      # Recreate Strata Vec to Pass In
      nstrat <- c(strataVecTrain, strataVecTest)
      
    }
    else
    {
      # Put all the train and test data in a single stratum
      distribution_obj$strata <- rep(1, nrow(gbm_data_obj$x))
      trainStrat <- c(num_train_rows, rep(NA, num_train_rows-1))
      testStrat <- c(num_test_rows, rep(NA, num_test_rows-1))
      nstrat <- c(trainStrat, testStrat)
    }
    
    # Sort response according to strata
    # i.order sets order of outputs
    if (attr(gbm_data_obj$y, "type") == "right")
    {
      sorted <- c(order(distribution_obj$strata[1:num_train_rows], -gbm_data_obj$y[1:num_train_rows, 1]),
                  order(distribution_obj$strata[(num_train_rows+1):nrow(gbm_data_obj$x)],
                        -gbm_data_obj$y[(num_train_rows+1):nrow(gbm_data_obj$x), 1])) 
      i_order <- c(order(distribution_obj$strata[1:num_train_rows], -gbm_data_obj$y[1:num_train_rows, 1]),
                   order(distribution_obj$strata[(num_train_rows+1):nrow(gbm_data_obj$x)],
                         -gbm_data_obj$y[(num_train_rows+1):nrow(gbm_data_obj$x), 1]) + num_train_rows)
    }
    else if (attr(gbm_data_obj$y, "type") == "counting") 
    {
      sorted <- cbind(c(order(distribution_obj$strata[1:num_train_rows], -gbm_data_obj$y[1:num_train_rows, 1]),
                        order(distribution_obj$strata[(num_train_rows+1):nrow(gbm_data_obj$x)],
                              -gbm_data_obj$y[(num_train_rows+1):nrow(gbm_data_obj$x), 1])),
                      c(order(distribution_obj$strata[1:num_train_rows], -gbm_data_obj$y[1:num_train_rows, 2]),
                        order(distribution_obj$strata[(num_train_rows+1):nrow(gbm_data_obj$x)],
                              -gbm_data_obj$y[(num_train_rows+1):nrow(gbm_data_obj$x), 2])))
      i_order <- c(order(distribution_obj$strata[1:num_train_rows], -gbm_data_obj$y[1:num_train_rows, 1]),
                   order(distribution_obj$strata[(num_train_rows+1):nrow(gbm_data_obj$x)], 
                         -gbm_data_obj$y[(num_train_rows+1):nrow(gbm_data_obj$x), 1]) + num_train_rows)
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