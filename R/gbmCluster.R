gbmCluster <- function(n, ncv, verbose){
    # If number of cores (n) not given, try to work it out from the number
    # that appear to be available and the number of CV folds.
    if (is.null(n)){
        n.cores <- detectCores()
        if (n.cores > 1){ n <- min(n.cores - 1, ncv) }
        else { n <- n.cores}
        if (verbose){
            cat("Detected", n.cores, "cores; will attempt to use", n, "\n")
        }
    }
    list(cluster=makeCluster(n), n.cores=n)
}
