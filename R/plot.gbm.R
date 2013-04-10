plot.gbm <- function(x,
                     i.var=1,
                     n.trees=x$n.trees,
                     continuous.resolution=100,
                     return.grid=FALSE,
                     type="link",
                     ...)
{
   if (!is.element(type, c("link", "response"))){
      stop( "type must be either 'link' or 'response'")
   }

   if(all(is.character(i.var)))
   {
      i <- match(i.var,x$var.names)
      if(any(is.na(i)))
      {
         stop("Plot variables not used in gbm model fit: ",i.var[is.na(i)])
      } else
      {
         i.var <- i
      }
   }

   if((min(i.var)<1) || (max(i.var)>length(x$var.names)))
   {
      warning("i.var must be between 1 and ",length(x$var.names))
   }
   if(n.trees > x$n.trees)
   {
      warning(paste("n.trees exceeds the number of trees in the model, ",x$n.trees,
                    ". Plotting using ",x$n.trees," trees.",sep=""))
      n.trees <- x$n.trees
   }

   if(length(i.var) > 3)
   {
      warning("gbm.int.plot creates up to 3-way interaction plots.\nplot.gbm will only return the plotting data structure.")
      return.grid = TRUE
   }

   # generate grid to evaluate gbm model
   grid.levels <- vector("list",length(i.var))
   for(i in 1:length(i.var))
   {
      # continuous
      if(is.numeric(x$var.levels[[i.var[i]]]))
      {
         grid.levels[[i]] <- seq(min(x$var.levels[[i.var[i]]]),
                                 max(x$var.levels[[i.var[i]]]),
                                 length=continuous.resolution)
      }
      # categorical or ordered
      else
      {
         grid.levels[[i]] <- as.numeric(factor(x$var.levels[[i.var[i]]],
                                               levels=x$var.levels[[i.var[i]]]))-1
      }
   }

   X <- expand.grid(grid.levels)
   names(X) <- paste("X",1:length(i.var),sep="")

   # Next if block for compatibility with objects created with 1.6
   if (is.null(x$num.classes)){
       x$num.classes <- 1
   }

   # evaluate at each data point
   y <- .Call("gbm_plot",
              X = as.double(data.matrix(X)),
              cRows = as.integer(nrow(X)),
              cCols = as.integer(ncol(X)),
              n.class = as.integer(x$num.classes),
              i.var = as.integer(i.var-1),
              n.trees = as.integer(n.trees) ,
              initF = as.double(x$initF),
              trees = x$trees,
              c.splits = x$c.splits,
              var.type = as.integer(x$var.type),
              PACKAGE = "gbm")

   if (x$distribution$name=="multinomial")
   {
      ## Put result into matrix form
      X$y <- matrix(y, ncol = x$num.classes)
      colnames(X$y) <- x$classes

      ## Use class probabilities
      if (type=="response"){
         X$y <- exp(X$y)
         X$y <- X$y / matrix(rowSums(X$y), ncol=ncol(X$y), nrow=nrow(X$y))
      }
   }
   else if(is.element(x$distribution$name, c("bernoulli", "pairwise")) && type=="response") {
      X$y <- 1/(1+exp(-y))
   }
   else if ((x$distribution$name=="poisson") && (type=="response")){
      X$y <- exp(y)
   }
   else if (type=="response"){
      warning("type 'response' only implemented for 'bernoulli', 'poisson', 'multinomial', and 'pairwise'. Ignoring" )
   }
   else { X$y <- y }

   # transform categorical variables back to factors
   f.factor <- rep(FALSE,length(i.var))
   for(i in 1:length(i.var))
   {
      if(!is.numeric(x$var.levels[[i.var[i]]]))
      {
         X[,i] <- factor(x$var.levels[[i.var[i]]][X[,i]+1],
                         levels=x$var.levels[[i.var[i]]])
         f.factor[i] <- TRUE
      }
   }

   if(return.grid)
   {
      names(X)[1:length(i.var)] <- x$var.names[i.var]
      return(X)
   }

   # create the plots
   if(length(i.var)==1)
   {
      if(!f.factor)
      {
         j <- order(X$X1)

         if (x$distribution$name == "multinomial") {
            if ( type == "response" ){
               ylabel <- "Predicted class probability"
            }
            else { ylabel <- paste("f(",x$var.names[i.var],")",sep="") }
            plot(range(X$X1), range(X$y), type = "n", xlab = x$var.names[i.var],
                 ylab = ylabel)

            for (ii in 1:x$num.classes){
               lines(X$X1,X$y[,ii],
                     xlab=x$var.names[i.var],
                     ylab=paste("f(",x$var.names[i.var],")",sep=""),
                     col = ii, ...)
            }
         }
         else if (is.element(x$distribution$name, c("bernoulli", "pairwise"))) {
            if ( type == "response" ){
               ylabel <- "Predicted probability"
            }
            else {
               ylabel <- paste("f(",x$var.names[i.var],")",sep="")
            }
            plot( X$X1, X$y , type = "l", xlab = x$var.names[i.var], ylab=ylabel )
         }
         else if ( x$distribution$name == "poisson" ){
            if (type == "response" ){
               ylabel <- "Predicted count"
            }
            else{
               ylabel <- paste("f(",x$var.names[i.var],")",sep="")
            }
            plot( X$X1, X$y , type = "l", xlab = x$var.names[i.var], ylab=ylabel )
         }
         else {
            plot(X$X1,X$y,
                 type="l",
                 xlab=x$var.names[i.var],
                 ylab=paste("f(",x$var.names[i.var],")",sep=""),...)
         }
      }
      else
      {
         if (x$distribution$name == "multinomial") {
            nX <- length(X$X1)
            dim.y <- dim(X$y)
            if (type == "response" ){
               ylabel <- "Predicted probability"
            }
            else{ ylabel <- paste("f(",x$var.names[i.var],")",sep="") }

            plot(c(0,nX), range(X$y), axes = FALSE, type = "n",
                 xlab = x$var.names[i.var], ylab = ylabel)
            axis(side = 1, labels = FALSE, at = 0:nX)
            axis(side = 2)

            mtext(as.character(X$X1), side = 1, at = 1:nX - 0.5)

            segments(x1 = rep(1:nX - 0.75, each = dim.y[2]), y1 = as.vector(t(X$y)),
                     x2 = rep(1:nX - 0.25, each = dim.y[2]), col = 1:dim.y[2])
         }
         else if (is.element(x$distribution$name, c("bernoulli", "pairwise")) && type == "response" ){
            ylabel <- "Predicted probability"
            plot( X$X1, X$y, type = "l", xlab=x$var.names[i.var], ylab=ylabel )
         }
         else if ( x$distribution$name == "poisson" & type == "response" ){
            ylabel <- "Predicted count"
            plot( X$X1, X$y, type = "l", xlab=x$var.names[i.var], ylab=ylabel )
         }
         else {
            plot(X$X1,X$y,
                 type="l",
                 xlab=x$var.names[i.var],
                 ylab=paste("f(",x$var.names[i.var],")",sep=""),...)
         }
      }
   }
   else if(length(i.var)==2)
   {
      if(!f.factor[1] && !f.factor[2])
      {
         if (x$distribution$name == "multinomial")
         {
            for (ii in 1:x$num.classes){
               X$temp <- X$y[, ii]
               print(levelplot(temp~X1*X2,data=X,
                               xlab=x$var.names[i.var[1]],
                               ylab=x$var.names[i.var[2]],...))
               title(paste("Class:", dimnames(X$y)[[2]][ii]))
            }
            X$temp <- NULL
         }
         else {
            print(levelplot(y~X1*X2,data=X,
                      xlab=x$var.names[i.var[1]],
                      ylab=x$var.names[i.var[2]],...))
         }
      }
      else if(f.factor[1] && !f.factor[2])
      {
         if (x$distribution$name == "multinomial")
         {
            for (ii in 1:x$num.classes){
               X$temp <- X$y[, ii]
               print( xyplot(temp~X2|X1,data=X,
                             xlab=x$var.names[i.var[2]],
                             ylab=paste("f(",x$var.names[i.var[1]],",",x$var.names[i.var[2]],")",sep=""),
                             type="l",
                             panel = panel.xyplot,
                             ...) )
               title(paste("Class:", dimnames(X$y)[[2]][ii]))
            }
            X$temp <- NULL
         }
         else {
            print(xyplot(y~X2|X1,data=X,
                   xlab=x$var.names[i.var[2]],
                   ylab=paste("f(",x$var.names[i.var[1]],",",x$var.names[i.var[2]],")",sep=""),
                   type="l",
                   panel = panel.xyplot,
                   ...))
         }
      }
      else if(!f.factor[1] && f.factor[2])
      {
         if (x$distribution$name == "multinomial")
         {
            for (ii in 1:x$num.classes){
               X$temp <- X$y[, ii]
               print( xyplot(temp~X1|X2,data=X,
                             xlab=x$var.names[i.var[1]],
                             ylab=paste("f(",x$var.names[i.var[1]],",",x$var.names[i.var[2]],")",sep=""),
                             type="l",
                             panel = panel.xyplot,
                             ...) )
               title(paste("Class:", dimnames(X$y)[[2]][ii]))
            }
            X$temp <- NULL
         }
         else {
            print(xyplot(y~X1|X2,data=X,
                   xlab=x$var.names[i.var[1]],
                   ylab=paste("f(",x$var.names[i.var[1]],",",x$var.names[i.var[2]],")",sep=""),
                   type="l",
                   panel = panel.xyplot,
                   ...))
         }
      }
      else
      {
         if (x$distribution$name == "multinomial")
         {
            for (ii in 1:x$num.classes){
               X$temp <- X$y[, ii]
               print( stripplot(X1~temp|X2,data=X,
                                xlab=x$var.names[i.var[2]],
                                ylab=paste("f(",x$var.names[i.var[1]],",",x$var.names[i.var[2]],")",sep=""),
                                ...) )
               title(paste("Class:", dimnames(X$y)[[2]][ii]))
            }
            X$temp <- NULL
         }
         else {
            print(stripplot(X1~y|X2,data=X,
                      xlab=x$var.names[i.var[2]],
                      ylab=paste("f(",x$var.names[i.var[1]],",",x$var.names[i.var[2]],")",sep=""),
                      ...))
         }
      }
   }
   else if(length(i.var)==3)
   {
      i <- order(f.factor)
      X.new <- X[,i]
      X.new$y <- X$y
      names(X.new) <- names(X)

      # 0 factor, 3 continuous
      if(sum(f.factor)==0)
      {
         X.new$X3 <- equal.count(X.new$X3)
         if (x$distribution$name == "multinomial")
         {
            for (ii in 1:x$num.classes){
               X.new$temp <- X.new$y[, ii]
               print( levelplot(temp~X1*X2|X3,data=X.new,
                                xlab=x$var.names[i.var[i[1]]],
                                ylab=x$var.names[i.var[i[2]]],...) )
               title(paste("Class:", dimnames(X.new$y)[[2]][ii]))
            }
            X.new$temp <- NULL
         }
         else {
            print(levelplot(y~X1*X2|X3,data=X.new,
                      xlab=x$var.names[i.var[i[1]]],
                      ylab=x$var.names[i.var[i[2]]],...))
         }
      }
      # 1 factor, 2 continuous
      else if(sum(f.factor)==1)
      {
         if (x$distribution$name == "multinomial")
         {
            for (ii in 1:x$num.classes){
               X.new$temp <- X.new$y[, ii]
               print( levelplot(temp~X1*X2|X3,data=X.new,
                                xlab=x$var.names[i.var[i[1]]],
                                ylab=x$var.names[i.var[i[2]]],...))
               title(paste("Class:", dimnames(X.new$y)[[2]][ii]) )
            }
            X.new$temp <- NULL
         }
         else {
            print(levelplot(y~X1*X2|X3,data=X.new,
                      xlab=x$var.names[i.var[i[1]]],
                      ylab=x$var.names[i.var[i[2]]],...))
         }
      }
      # 2 factors, 1 continuous
      else if(sum(f.factor)==2)
      {
         if (x$distribution$name == "multinomial")
         {
            for (ii in 1:x$num.classes){
               X.new$temp <- X.new$y[, ii]
               print( xyplot(temp~X1|X2*X3,data=X.new,
                             type="l",
                             xlab=x$var.names[i.var[i[1]]],
                             ylab=paste("f(",paste(x$var.names[i.var[1:3]],collapse=","),")",sep=""),
                             panel = panel.xyplot,
                             ...) )
               title(paste("Class:", dimnames(X.new$y)[[2]][ii]) )
            }
            X.new$temp <- NULL
         }
         else {
            print(xyplot(y~X1|X2*X3,data=X.new,
                   type="l",
                   xlab=x$var.names[i.var[i[1]]],
                   ylab=paste("f(",paste(x$var.names[i.var[1:3]],collapse=","),")",sep=""),
                   panel = panel.xyplot,
                   ...))
         }
      }
      # 3 factors, 0 continuous
      else if(sum(f.factor)==3)
      {
         if (x$distribution$name == "multinomial")
         {
            for (ii in 1:x$num.classes){
               X.new$temp <- X.new$y[, ii]
               print( stripplot(X1~temp|X2*X3,data=X.new,
                                xlab=x$var.names[i.var[i[1]]],
                                ylab=paste("f(",paste(x$var.names[i.var[1:3]],collapse=","),")",sep=""),
                                ...) )
               title(paste("Class:", dimnames(X.new$y)[[2]][ii]) )
            }
            X.new$temp <- NULL
         }
         else {
            print(stripplot(X1~y|X2*X3,data=X.new,
                      xlab=x$var.names[i.var[i[1]]],
                      ylab=paste("f(",paste(x$var.names[i.var[1:3]],collapse=","),")",sep=""),
                      ...))
         }
      }
   }
}
