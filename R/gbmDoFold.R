gbmDoFold <- function(X,
         i.train, x, y, offset, distribution, w, var.monotone, n.trees,
         interaction.depth, n.minobsinnode, shrinkage, bag.fraction, mFeatures,
         cv.group, var.names, response.name, group, s, lVerbose, keep.data, nTrain, tied.times.method, prior.node.coeff.var,
         strata){
    # Do specified cross-validation fold - a self-contained function for
    # passing to individual cores.

    library(gbm, quietly=TRUE)
   
    # Handle the final model case separately
    if (X == 0){
      if (lVerbose) message("Fitting Final Model \n")
      res <- gbm.fit(x,y,
                       offset = offset,
                       distribution = distribution,
                       w = w,
                       var.monotone = var.monotone,
                       n.trees = n.trees,
                       interaction.depth = interaction.depth,
                       n.minobsinnode = n.minobsinnode,
                       shrinkage = shrinkage,
                       bag.fraction = bag.fraction,
                       nTrain = nTrain,
                       mFeatures = mFeatures,
                       keep.data = keep.data,
                       verbose = lVerbose,
                       var.names = var.names,
                       response.name = response.name,
                       group = group,
                      misc = tied.times.method,
                     prior.node.coeff.var = prior.node.coeff.var,
                     strata = strata)
    } else {
      if (lVerbose) message("CV:", X, "\n")
      set.seed(s[[X]])
      i <- order(cv.group == X)
      x <- x[i.train,,drop=FALSE][i,,drop=FALSE]
      y <- y[i.train][i]
      offset <- offset[i.train][i]
      nTrain <- length(which(cv.group != X))
      group <- group[i.train][i]
      
      res <- gbm.fit(x, y,
                     offset=offset, distribution=distribution,
                     w=w, var.monotone=var.monotone, n.trees=n.trees,
                     interaction.depth=interaction.depth,
                     n.minobsinnode=n.minobsinnode,
                     shrinkage=shrinkage,
                     bag.fraction=bag.fraction,
                     nTrain=nTrain, mFeatures=mFeatures, keep.data=FALSE,
                     verbose=FALSE, response.name=response.name,
                     group=group, misc=tied.times.method, prior.node.coeff.var=prior.node.coeff.var,
                     strata = strata)
  }
  res
}
