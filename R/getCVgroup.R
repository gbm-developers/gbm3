getCVgroup <-
    # Construct cross-validation groups depending on the type of model to be fit
function(distribution, class.stratify.cv, y, i.train, cv.folds, group){

    if (distribution$name %in% c( "bernoulli", "multinomial" ) & class.stratify.cv ){
        nc <- table(y[i.train]) # Number in each class
        uc <- names(nc)
        if (min(nc) < cv.folds){
            stop( paste("The smallest class has only", min(nc), "objects in the training set. Can't do", cv.folds, "fold cross-validation."))
        }
        cv.group <- vector(length = length(i.train))
        for (i in 1:length(uc)){
           cv.group[y[i.train] == uc[i]] <- sample(rep(1:cv.folds , length = nc[i]))
        }
    } # Close if
    else if (distribution$name == "pairwise") {
         # Split into CV folds at group boundaries
         s <- sample(rep(1:cv.folds, length=nlevels(group)))
         cv.group <- s[as.integer(group[i.train])]
      }
      else {
         cv.group <- sample(rep(1:cv.folds, length=length(i.train)))
      }
      cv.group
}
