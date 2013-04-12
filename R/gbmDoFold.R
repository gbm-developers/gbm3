gbmDoFold <-
    # Do specified cross-validation fold - a self-contained function for
    # passing to individual cores.
function(X,
         i.train, x, y, offset, distribution, w, var.monotone, n.trees,
         interaction.depth, n.minobsinnode, shrinkage, bag.fraction,
         cv.group, var.names, response.name, group, s){
    library(gbm, quietly=TRUE)
    cat("CV:", X, "\n")

    set.seed(s[[X]])

    i <- order(cv.group == X)
    x <- x[i.train,,drop=TRUE][i,,drop=FALSE]
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
                   nTrain=nTrain, keep.data=FALSE,
                   verbose=FALSE, response.name=response.name,
                   group=group)
    res
}
