gbmDoFold <- function(X, i.train, x, y, offset, distribution, w,
                      var.monotone, n.trees, interaction.depth,
                      n.minobsinnode, shrinkage, bag.fraction,
                      mFeatures, cv.group, var.names, response.name,
                      group, s, lVerbose, keep.data, nTrain,
                      tied.times.method, prior.node.coeff.var, strata,
                      patient.id, par.details){
    # Do specified cross-validation fold - a self-contained function for
    # passing to individual cores.

    library(gbm, quietly=TRUE)
   
    # Handle the final model case separately
    if (X == 0){
      if (lVerbose) message("Fitting Final Model \n")
      res <- gbm.fit(x,y, offset = offset,
                     distribution = distribution, w = w,
                     var.monotone = var.monotone, n.trees = n.trees,
                     interaction.depth = interaction.depth,
                     n.minobsinnode = n.minobsinnode,
                     shrinkage = shrinkage,
                     bag.fraction = bag.fraction, nTrain = nTrain,
                     mFeatures = mFeatures, keep.data = keep.data,
                     verbose = lVerbose, var.names = var.names,
                     response.name = response.name, group = group,
                     misc = tied.times.method,
                     prior.node.coeff.var = prior.node.coeff.var,
                     strata = strata, patient.id = patient.id,
                     par.details = par.details)
    } else {
      if (lVerbose) message("CV:", X, "\n")
      set.seed(s[[X]])
      
      # Patients in training set
      patients_in_training_set <- patient.id %in% i.train
      patient.id <- patient.id[patients_in_training_set]
      
      # Patients in cv.group
      patient_id_in_cv_group <- patient.id %in% i.train[(cv.group == X)]
            
      # Pulls out patients in cv.group
      i <- order(patient_id_in_cv_group)
      
      # Get the correct cv group 
      x <- x[patients_in_training_set,,drop=FALSE][i,,drop=FALSE]
      y <- y[patients_in_training_set][i]
      offset <- offset[patients_in_training_set][i]
      nTrain <- length(which(cv.group != X))
      group <- group[patients_in_training_set][i]
      
      
      if(distribution$name=="coxph") {
        strata <- strata[patients_in_training_set][i]
      }
      
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
                     strata = strata, patient.id = patient.id, par.details=par.details)
  }
  res
}
