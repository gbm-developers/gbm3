gbmCluster <- function(n){
    if (!is.element('package:parallel', search())){
        cat("Attaching 'parallel' package. ")
        library(parallel)
    }
    if (is.null(n)){
        n <- detectCores()
        cat("Detected", n, "cores.\n")
    }
    list(cluster=makeCluster(n), n.cores=n)
}
