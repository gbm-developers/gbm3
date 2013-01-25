getStratify <- function(strat, d){
    if (is.null(strat)){
        if (d$name == "multinomial" ){ strat <- TRUE }
        else { strat <- FALSE }
    }
    else {
        if (!is.element(d$name, c( "bernoulli", "multinomial"))){
           warning("You can only use class.stratify.cv when distribution is bernoulli or multinomial. Ignored.")
           strat <- FALSE
        }
    } # Close else
    strat
}
