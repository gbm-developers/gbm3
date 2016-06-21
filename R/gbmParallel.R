##' Control parallelization options
##'
##' \code{gbm} uses openmp to parallelize its core algorithm, and the
##' details are controlled by this object.  As guidance, set
##' \code{n.threads} to the number of cores on your computer, and set
##' \code{arrayChunkSize} to a reasonable - not necessarily small -
##' size.
##' 
##' @param n.threads the number of threads to use (a positive
##'     integer).  The number of cores on your computer is a
##'     reasonable default.
##' @param arrayChunkSize the size of chunks to use in array scans.
##'     Values that are too small result in a great deal of overhead;
##'     The default of 1024 appears reasonable, but do experiment.
##' @return an object of type \code{gbmParallel}
##' @export
gbmParallel <- function(n.threads=1, arrayChunkSize=1024) {
    res <- list(n.threads=n.threads,
                arrayChunkSize=arrayChunkSize)
    class(res) <- "gbmParallel"
    res
}

##' @export
print.gbmParallel <- function(x, ...) {
    cat("GBM parallelization\n\n",
        "number of threads: ", x$n.threads, "\n",
        "array chunk size : ", x$arrayChunkSize, "\n",
        sep="")
    invisible(x)
}

