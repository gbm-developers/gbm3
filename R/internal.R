## miscellaneous small internal helper functions.

## these are not exported and not formally documented in the manual

## Warn if a variable does not vary
##
## @param x a numeric variable to check
## @param ind an index to use in a potential warning
## @param name a name to include in a potential warning
## @return NULL, invisibly

warnNoVariation <- function(x, ind, name) {
    ## suppress warnings in here
    ## because min and max warn if the variable is completely NA
    suppressWarnings(variation <- range(x, na.rm=TRUE))

    ## I really mean ">=" here, which catches the all NA case
    ## and the standard case
    if (variation[[1]] >= variation[[2]]) {
        warning("variable ", ind, ": ", name, " has no variation.")
    }

    invisible(NULL)
}
