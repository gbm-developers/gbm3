#' Training Groups
#' 
#' Determine the number of training groups - update
#' weights and input data ordering
#' 
#' @usage  training_groups_reorder(x, y, weights, feature_names, training_params, distribution_obj)
#' 
#' @param x 
#' 
#' @param y
#' 
#' @param weights
#' 
#' @param feature_names
#' 
#' @param train_params
#' 
#' @param distribution_obj
#' 
#' @return a list containing the reordered input data, updated weights, training parameters
#' and distribution objects
#' 

training_groups_reorder <- function(x, y, weights, feature_names, distribution_obj) {
  if(distribution_obj$reorder) {
   UseMethod("trainig_groups_reorder", distribution_obj) 
  }
}

training_groups_reorder.CoxPHGBMDist <- function(x, y, weights, feature_names, train_params,
                                                 distribution_obj) {
  
}

training_groups_reorder.PairwiseGBMDist <- function(x, y, weights, feature_names, train_params, 
                                                    distribution_obj) {
  # Set up output
  output <- list()
  
  if (is.null(distribution.group))
  {
    stop("For pairwise regression, the distribution parameter must be a list with a parameter 'group' for the list of the column names indicating groups, for example list(name=\"pairwise\",group=c(\"date\",\"session\",\"category\",\"keywords\")).")
  }
  
  # Check if group names are valid
  i <- match(distribution_obj$group, feature_names)
  if (any(is.na(i)))
  {
    stop("Group column does not occur in data: ", distribution_obj$group[is.na(i)])
  }
  
  # Construct group index
  distribution_obj$group <- factor(do.call(paste, c(data[,distribution_obj$group, drop=FALSE], sep=":")))
  
  # Check that weights are constant across groups
  if ((!missing(weights)) && (!is.null(weights)))
  {
    weights.min <- tapply(weights, INDEX=group, FUN=min)
    weights.max <- tapply(weights, INDEX=group, FUN=max)
    
    if (any(weights.min != w.max))
    {
      stop("For distribution 'pairwise', all instances for the same group must have the same weight")
    }
    
    # Normalize across groups
    weights <- weights * length(w.min) / sum(w.min)
  }
  
  # Shuffle groups, to remove bias when splitting into train/test set and/or CV folds
  perm.levels  <- levels(distribution_obj$group)[sample(1:nlevels(distribution_obj$group))]
  distribution_obj$group        <- factor(distribution_obj$group, levels=perm.levels)
  
  # The C function expects instances to be sorted by group and descending by target
  ord.group    <- order(distribution_obj$group, -y)
  distribution_obj$group  <- group[ord.group]
  y            <- y[ord.group]
  x            <- x[ord.group,,drop=FALSE]
  w            <- w[ord.group]
  
  # Split into train and validation set, at group boundary
  train_params$num_train <- max(1, round(train.fraction * nlevels(distribution_obj$group)))
}