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
