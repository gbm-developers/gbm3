## Accessor methods for gbm fit objects

##' Extract trees from GBMFit objects
##'
##' @param obj a GBMFit object
##' @return a vector containing the trees generated in the fit
##' @export
trees <- function(obj) {
    check_if_gbm_fit(obj)
    obj$trees
}

##' Extract errors from GBMFit objects
##'
##' @param obj a GBMFit object
##' @param which which error measure to extract
##' @return a vector giving the error by tree
##' @export
iteration_error <- function(obj, which=c('train', 'valid', 'cv')) {
    check_if_gbm_fit(obj)
    switch(match.arg(which),
           train=obj$train.error,
           valid=obj$valid.error,
           cv=obj$cv_error,
           stop("Unknown error measure"))
}

##' What is the distribution name used here?
##'
##' @param obj the object for which the distribution name is needed
##' @param ... other parameters
##' @return a string identifying the distribution
##' @export
distribution_name <- function(obj, ...) {
    UseMethod("distribution_name")
}

##' @export
distribution_name.default <- function(obj, ...) {
    stop("I don't know how to handle this")
}

##' @export
distribution_name.GBMFit <- function(obj, ...) {
    distribution_name(obj$distribution)
}

##' @export
distribution_name.GBMDist <- function(obj, ...) {
    obj$name
}
