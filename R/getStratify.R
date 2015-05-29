getStratify <- function(strat, d){
    if (is.null(strat)){
        strat <- FALSE
    }
    else {
        if (!(d$name %in% c("bernoulli"))){
            warning("You can only use class.stratify.cv when distribution is bernoulli. Ignored.")
        }
        strat <- FALSE
    }
    strat
}
