#' @describeIn gbm add new trees to a gbm object
#' @export
gbm.more <- function(object,
                     n.new.trees = 100,
                     data = NULL,
                     weights = NULL,
                     offset = NULL,
                     verbose = NULL,
                     par.details = getOption("gbm.parallel"))
{
   return(object)
}
