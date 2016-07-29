#' Performance for pairwise
#' 
#' Additional performance using appropriate metric for pairwise distribution fit.
#' 
#' @usage perf_pairwise(y, f, group, metric='ndcg', w=NULL, max_rank=0)
#' 
#' @param y responses used for fit.
#' 
#' @param f the predicted responses.
#' 
#' @param metric What type of performance measure to compute in \code{perf_pairwise}.
#'   Can take values "ir_measure_conc", "ir_measure_mrr", "ir_measure_map" or
#'   "ir_measure_ndgc".
#' @param w is the weights of each observation.
#' 
#' @param max_rank the maximum rank considered in the groups measure. 

#' @return \code{gbm_perf} returns the estimated optimal number of iterations.
#' The method of computation depends on the \code{method} argument.
#' 
#' @author Greg Ridgeway \email{gregridgeway@@gmail.com}
#' 
#' @seealso \code{\link{gbm2}}
#' @keywords nonlinear survival nonparametric tree
#' @export

perf_pairwise <- function(y, f, group, metric="ndcg", w=NULL, max_rank=0){
  func.name <- switch(metric,
                      conc = "ir_measure_conc",
                      mrr  = "ir_measure_mrr",
                      map  = "ir_measure_map",
                      ndcg = "ir_measure_ndcg",
                      stop(paste("Metric",metric,"is not supported"))
  )
  
  # Optimization: for binary targets,
  # AUC is equivalent but faster than CONC
  if (metric == "conc" && all(is.element(y, 0:1))) {
    func.name <- "ir_measure_auc"
  }
  
  # Max rank = 0 means no cut off
  if (max_rank <= 0) {
    max_rank <- length(y)+1
  }
  
  # Random tie breaking in case of duplicate scores.
  # (Without tie breaking, we would overestimate if instances are
  # sorted descending on target)
  f <- f + 1E-10 * runif(length(f), min=-0.5, max=0.5)
  
  measure.by.group <- as.matrix(by(list(y, f), INDICES=group, FUN=get(func.name), max_rank=max_rank))
  
  # Exclude groups with single result or only negative or positive instances
  idx <- which((!is.null(measure.by.group)) & measure.by.group >= 0)
  
  if (is.null(w)) {
    return (mean(measure.by.group[idx]))
  } else {
    # Assumption: weights are constant per group
    w.by.group <- tapply(w, group, mean)
    return (weighted.mean(measure.by.group[idx], w=w.by.group[idx]))
  }
}


#### Helper Functions ####
ir_measure_conc <- function(y.f, max_rank=0) {
  # Note: max_rank is meaningless for CONC
  
  y           <- y.f[[1]]
  f           <- y.f[[2]]
  
  tab         <- table(y)
  csum        <- cumsum(tab)
  total.pairs <- sum(tab * (csum - tab))
  
  if (total.pairs == 0) {
    return (-1.0)
  } else {
    return (gbm_conc(y[order(-f)]) / total.pairs)
  }
}

ir_measure_auc <- function(y.f, max_rank=0){
  # Note: max_rank is meaningless for AUC
  y       <- y.f[[1]]
  f       <- y.f[[2]]
  num.pos <- sum(y>0)
  
  if (length(f) <= 1 || num.pos == 0 || num.pos == length(f))
  {
    return (-1.0)
  }
  else
  {
    return (gbm_roc_area(obs=y, pred=f))
  }
}

ir_measure_mrr <- function(y.f, max_rank) {
  y       <- y.f[[1]]
  f       <- y.f[[2]]
  num.pos <- sum(y>0)
  
  if (length(f) <= 1 || num.pos == 0 || num.pos == length(f))
  {
    return (-1.0)
  }
  
  ord         <- order(f, decreasing=TRUE)
  min.idx.pos <- min(which(y[ord]>0))
  
  if (min.idx.pos <= max_rank)
  {
    return (1.0 / min.idx.pos)
  }
  else
  {
    return (0.0)
  }
}

ir_measure_map <- function(y.f, max_rank=0) {
  # Note: max_rank is meaningless for MAP
  
  y         <- y.f[[1]]
  f         <- y.f[[2]]
  ord       <- order(f, decreasing=TRUE)
  idx.pos   <- which(y[ord]>0)
  num.pos   <- length(idx.pos)
  
  if (length(f) <= 1 || num.pos == 0 || num.pos == length(f))
  {
    return (-1.0)
  }
  
  # Above and including the rank of the i-th positive result,
  # there are exactly i positives and rank(i) total results
  return (sum((1:length(idx.pos))/idx.pos) / num.pos)
}

ir_measure_ndcg <- function(y.f, max_rank) {
  y         <- y.f[[1]]
  f         <- y.f[[2]]
  
  if (length(f) <= 1 || all(diff(y)==0)) return (-1.0)
  
  num.items <- min(length(f), max_rank)
  ord       <- order(f, decreasing=TRUE)
  
  dcg       <- sum(y[ord][1:num.items] / log2(2:(num.items+1)))
  
  # The best possible DCG: order by target
  ord.max   <- order(y, decreasing=TRUE)
  dcg.max   <- sum(y[ord.max][1:num.items] / log2(2:(num.items+1)))
  
  # Normalize
  return (dcg / dcg.max)
}