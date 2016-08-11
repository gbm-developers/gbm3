# Empty Distribution
#  
# Creates an empty distribution object
# 
# @usage empty_distribution(name)
# 
# @param name A string specifying which empty distribution object to build.
# 
# @author James Hickey
#
# @return An empty distribution class - only the name defined.
#

empty_distribution <- function(name) {
  return(structure(list(name=name, reorder=FALSE),
                   class=c(paste0(name, "GBMDist"), "GBMDist")))
}