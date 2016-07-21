#' Empty Distribution
#'  
#' Creates an empty distribution object
#' 
#' @usage empty_distribution(name)
#' 
#' @param name A string specifying which empty distribution object to build.
#' 
#' @return An empty distribution class - only the name defined.
#' 
#' @export

empty_distribution <- function(name) {
  return(structure(list(name=name, reorder=FALSE),
                   class=c(paste0(name, "GBMDist"), "GBMDist")))
}