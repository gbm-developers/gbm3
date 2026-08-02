# Delete observations that are missing a response or missing all of the 
#   predictors, anything else is kept
nagbm <- function(x){
    Terms <- attr(x, 'terms')
    if(!is.null(Terms)) yvar <- attr(Terms, "response") else yvar <- 0L
    if (yvar==0L) {
	remove <- apply(is.na(x), 1, all)
	}
    else {
	xmiss <- is.na(x[-yvar])
	ymiss <- is.na(x[[yvar]])
	if (is.matrix(ymiss))
	    remove <- (apply(xmiss, 1, all) | apply(ymiss, 1, any)) 
	else
	    remove <- (apply(xmiss, 1, all) | ymiss)
	}
    if (!any(remove)) x
    else {
	temp <- seq(remove)[remove]   # list of dropped rows
	names(temp) <- row.names(x)[remove]
	#the methods for this group are all the same as for na.omit
	class(temp) <- c("nagbm", "omit")
	structure(x[!remove,,drop=FALSE], na.action=temp)
	}
    }
