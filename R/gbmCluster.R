gbmCluster <- function(n){
    # If number of cores (n) not given, try to work it out from the number
    # that appear to be available and the number of CV folds.
    if (is.null(n)){
      n <- detectCores()
    }
    if (n == 1){ NULL }
    else { makeCluster(n) }
}
