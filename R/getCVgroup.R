# Construct cross-validation groups 
# depending on the type of model to be fit

getCVgroup <- function(
  distribution, class.stratify.cv, y
  , i.train, cv.folds, group, fold.id){
  
  if(distribution == "bernoulli" & class.stratify.cv){
    
    # Number in each class
    Ones <- tabulate(y[i.train])
    Zeros <- length(y[i.train])-Ones
    
    smallGroup <- min(c(Ones, Zeros))
    
    if(smallGroup < cv.folds){
      stop(
        paste("The smallest class has only"
          ,smallGroup
          ,"objects in the training set. Can't do"
          ,cv.folds, "fold cross-validation.")
        )
      }
    
    cv.group <- c(
      sample(rep(1:cv.folds, length = Zeros))
    , sample(rep(1:cv.folds, length = Ones))
    )
    
  } else if (distribution == "pairwise") {
    
    # Split into CV folds at group boundaries
    s <- sample(rep(1:cv.folds, length=nlevels(group)))
    cv.group <- s[as.integer(group[i.train])]
         
  } else if (!is.null(fold.id)) {
    
    cv.group <- fold.id
    
  } else {
    
    cv.group <- sample(rep(1:cv.folds, length=length(i.train)))

  }
  
  cv.group
  
}
