##' Control parallelization options
##'
##' \code{gbm} uses openmp to parallelize its core algorithm, and the
##' details are controlled by this object.  As guidance, set
##' \code{n.threads} to the number of cores on your computer, and set
##' \code{arrayChunkSize} to a reasonable - not necessarily small -
##' size.
##' 
##' @param num_threads the number of threads to use (a positive
##'     integer).  The number of cores on your computer is a
##'     reasonable default.
##' @param array_chunk_size the size of chunks to use in array scans.
##'     Values that are too small result in a great deal of overhead;
##'     The default of 1024 appears reasonable, but do experiment.
##' @return an object of type \code{gbmParallel}
##' @export
gbmParallel <- function(num_threads=1, array_chunk_size=1024) {
    res <- list(num_threads=num_threads,
                array_chunk_size=array_chunk_size)
    class(res) <- "gbmParallel"
    res
}

##' @export
print.gbmParallel <- function(x, ...) {
    cat("GBM parallelization\n\n",
        "number of threads: ", x$num_threads, "\n",
        "array chunk size : ", x$array_chunk_size, "\n",
        sep="")
    invisible(x)
}

