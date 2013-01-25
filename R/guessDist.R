guessDist <- function(y){
    # If distribution is not given, try to guess it
    if (length(unique(y)) == 2){ d <- "bernoulli" }
    else if (class(y) == "Surv" ){ d <- "coxph" }
    else if (is.factor(y)){ d <- "multinomial" }
    else{ d <- "gaussian" }
    cat(paste("Distribution not specified, assuming", d, "...\n"))
    list(name=d)
}
