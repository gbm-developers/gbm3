# Compute Friedman's H statistic for interaction effects


#' Estimate the strength of interaction effects
#' 
#' Computes Friedman's H-statistic to assess the strength of variable
#' interactions.
#' 
#' \code{interact.gbm} computes Friedman's H-statistic to assess the relative
#' strength of interaction effects in non-linear models. H is on the scale of
#' [0-1] with higher values indicating larger interaction effects. To connect
#' to a more familiar measure, if \eqn{x_1} and \eqn{x_2} are uncorrelated
#' covariates with mean 0 and variance 1 and the model is of the form
#' \deqn{y=\beta_0+\beta_1x_1+\beta_2x_2+\beta_3x_3} then
#' \deqn{H=\frac{\beta_3}{\sqrt{\beta_1^2+\beta_2^2+\beta_3^2}}}
#' 
#' Note that if the main effects are weak, the estimated H will be unstable.
#' For example, if (in the case of a two-way interaction) neither main effect
#' is in the selected model (relative influence is zero), the result will be
#' 0/0. Also, with weak main effects, rounding errors can result in values of H
#' > 1 which are not possible.
#' 
#' @param x a \code{\link{gbm.object}} fitted using a call to \code{\link{gbm}}
#' @param data the dataset used to construct \code{x}. If the original dataset
#' is large, a random subsample may be used to accelerate the computation in
#' \code{interact.gbm}
#' @param i.var a vector of indices or the names of the variables for compute
#' the interaction effect. If using indices, the variables are indexed in the
#' same order that they appear in the initial \code{gbm} formula.
#' @param n.trees the number of trees used to generate the plot. Only the first
#' \code{n.trees} trees will be used
#' @return Returns the value of \eqn{H}.
#' @author Greg Ridgeway \email{gregridgeway@@gmail.com}
#' @seealso \code{\link{gbm}}, \code{\link{gbm.object}}
#' @references J.H. Friedman and B.E. Popescu (2005). \dQuote{Predictive
#' Learning via Rule Ensembles.} Section 8.1
#' @keywords methods
interact.gbm <- function(x, data, i.var = 1, n.trees = x$n.trees){
   ###############################################################
   # Do sanity checks on the call
    if (x$interaction.depth < length(i.var)){
       stop("interaction.depth too low in model call")
   }

   if (all(is.character(i.var))){
      i <- match(i.var, x$var.names)
      if (any(is.na(i))) {
         stop("Variables given are not used in gbm model fit: ", i.var[is.na(i)])
      }
      else {
         i.var <- i
      }
   }
   if ((min(i.var) < 1) || (max(i.var) > length(x$var.names))) {
      warning("i.var must be between 1 and ", length(x$var.names))
   }
   if (n.trees > x$n.trees) {
      warning(paste("n.trees exceeds the number of trees in the model, ",
                    x$n.trees,". Using ", x$n.trees, " trees.", sep = ""))
      n.trees <- x$n.trees
   }
   # End of sanity checks
   ###############################################################

   unique.tab <- function(z,i.var) {
      a <- unique(z[,i.var,drop=FALSE])
      a$n <- table(factor(apply(z[,i.var,drop=FALSE],1,paste,collapse="\r"),
                          levels=apply(a,1,paste,collapse="\r")))
      return(a)
   }

   # convert factors
   for(j in i.var) {
      if(is.factor(data[,x$var.names[j]]))
      data[,x$var.names[j]] <-
      as.numeric(data[,x$var.names[j]])-1
   }

   # generate a list with all combinations of variables
   a <- apply(expand.grid(rep(list(c(FALSE,TRUE)), length(i.var)))[-1,],1,
              function(x) as.numeric(which(x)))
   FF <- vector("list",length(a))
   for(j in 1:length(a)) {
      FF[[j]]$Z <- data.frame(unique.tab(data, x$var.names[i.var[a[[j]]]]))
      FF[[j]]$n <- as.numeric(FF[[j]]$Z$n)
      FF[[j]]$Z$n <- NULL
      FF[[j]]$f <- .Call("gbm_plot",
                         X = data.matrix(FF[[j]]$Z),
                         n.class = as.integer(x$num.classes),
                         i.var = as.integer(i.var[a[[j]]] - 1),
                         n.trees = as.integer(n.trees),
                         initF = as.double(x$initF),
                         trees = x$trees,
                         c.splits = x$c.splits,
                         var.type = as.integer(x$var.type),
                         PACKAGE = "gbm")
      # FF[[jj]]$Z is the data, f is the predictions, n is the number of levels for factors

      # Need to restructure f to deal with multinomial case
      FF[[j]]$f <- matrix(FF[[j]]$f, ncol=x$num.classes, byrow=FALSE)

      # center the values
      FF[[j]]$f <- apply(FF[[j]]$f, 2, function(x, w){
                                          x - weighted.mean(x, w, na.rm=TRUE)
                                       }, w=FF[[j]]$n)

      # precompute the sign of these terms to appear in H
      FF[[j]]$sign <- ifelse(length(a[[j]]) %% 2 == length(i.var) %% 2, 1, -1)
   }

   H <- FF[[length(a)]]$f

   for(j in 1:(length(a)-1)){
      i1 <- apply(FF[[length(a)]]$Z[,a[[j]], drop=FALSE], 1, paste, collapse="\r")
      i2 <- apply(FF[[j]]$Z,1,paste,collapse="\r")
      i <- match(i1, i2)

      H <- H + with(FF[[j]], sign*f[i,])
   }

   # Compute H
   w <- matrix(FF[[length(a)]]$n, ncol=1)
   f <- matrix(FF[[length(a)]]$f^2, ncol=x$num.classes, byrow=FALSE)

   top <- apply(H^2, 2, weighted.mean, w = w, na.rm = TRUE)
   btm <- apply(f, 2, weighted.mean, w = w, na.rm = TRUE)
   H <- top / btm

   if (x$distribution$name=="multinomial"){
      names(H) <- x$classes
   }

   # If H > 1, rounding and tiny main effects have messed things up
   H[H > 1] <- NaN

   return(sqrt(H))
}
