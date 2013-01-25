# Compute Friedman's H statistic for interaction effects
interact.gbm <- function(x, data, i.var = 1, n.trees = x$n.trees)
{
   if (all(is.character(i.var)))
   {
      i <- match(i.var, x$var.names)
      if (any(is.na(i))) {
         stop("Variables given are not used in gbm model fit: ",
              i.var[is.na(i)])
      }
      else
      {
         i.var <- i
      }
   }
   if ((min(i.var) < 1) || (max(i.var) > length(x$var.names)))
   {
      warning("i.var must be between 1 and ", length(x$var.names))
   }
   if (n.trees > x$n.trees)
   {
      warning(paste("n.trees exceeds the number of trees in the model, ",
                    x$n.trees,". Using ", x$n.trees, " trees.", sep = ""))
      n.trees <- x$n.trees
   }

   unique.tab <- function(z,i.var)
   {
      a <- unique(z[,i.var,drop=FALSE])
      a$n <- table(factor(apply(z[,i.var,drop=FALSE],1,paste,collapse="\r"),
                          levels=apply(a,1,paste,collapse="\r")))
      return(a)
   }

   # convert factors
   for(j in i.var)
   {
      if(is.factor(data[,x$var.names[j]]))
      data[,x$var.names[j]] <-
      as.numeric(data[,x$var.names[j]])-1
   }

   # generate a list with all combinations of variables
   a <- apply(expand.grid(rep(list(c(FALSE,TRUE)), length(i.var)))[-1,],1,
              function(x) as.numeric(which(x)))
   FF <- vector("list",length(a))
   for(j in 1:length(a))
   {
      FF[[j]]$Z <- data.frame(unique.tab(data, x$var.names[i.var[a[[j]]]]))
      FF[[j]]$n <- as.numeric(FF[[j]]$Z$n)
      FF[[j]]$Z$n <- NULL
      FF[[j]]$f <- .Call("gbm_plot",
                         X = as.double(data.matrix(FF[[j]]$Z)),
                         cRows = as.integer(nrow(FF[[j]]$Z)),
                         cCols = as.integer(ncol(FF[[j]]$Z)),
                         n.class = as.integer(x$num.classes),
                         i.var = as.integer(i.var[a[[j]]] - 1),
                         n.trees = as.integer(n.trees),
                         initF = as.double(x$initF),
                         trees = x$trees,
                         c.splits = x$c.splits,
                         var.type = as.integer(x$var.type),
                         PACKAGE = "gbm")
      # center the values
      FF[[j]]$f <- with(FF[[j]], f - weighted.mean(f,n))
      # precompute the sign of these terms to appear in H
      FF[[j]]$sign <- ifelse(length(a[[j]]) %% 2 == length(i.var) %% 2, 1, -1)
   }

   H <- FF[[length(a)]]$f
   for(j in 1:(length(a)-1))
   {
      i <- match(apply(FF[[length(a)]]$Z[,a[[j]],drop=FALSE],1,paste,collapse="\r"),
                 apply(FF[[j]]$Z,1,paste,collapse="\r"))
      H <- H + with(FF[[j]], sign*f[i])
   }
   if (is.null(dim(H)))
   {
      H <- weighted.mean(H^2, FF[[length(a)]]$n)/
              weighted.mean((FF[[length(a)]]$f)^2,FF[[length(a)]]$n)
   }
   else {
      H <- apply(H^2, 2, weighted.mean, w = FF[[length(a)]]$n, na.rm = TRUE)/
              apply((FF[[length(a)]]$f)^2, 2, weighted.mean,
                 w = FF[[length(a)]]$n, na.rm = TRUE)
   }
   if (x$distribution$name=="multinomial"){
      names(H) <- x$classes
   }
   return(sqrt(H))
}
