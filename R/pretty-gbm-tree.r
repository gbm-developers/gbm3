#' Print gbm tree components
#' 
#' \code{GBMFit} stores the collection of trees used to construct the model in a
#' compact matrix structure. This function extracts the information from a
#' single tree and displays it in a slightly more readable form. This function
#' is mostly for debugging purposes and to satisfy some users' curiosity.
#' 
#' @param gbm_fit_obj a \code{GBMFit} object initially fit using
#' \code{\link{gbm2}}.
#' 
#' @param tree_index the index of the tree component to extract from \code{gbm_fit_obj}
#' and display.
#' 
#' @return \code{pretty.GBMFit} returns a data frame. Each row corresponds to
#' a node in the tree. Columns indicate \item{SplitVar}{index of which variable
#' is used to split. -1 indicates a terminal node.} \item{SplitCodePred}{if the
#' split variable is continuous then this component is the split point. If the
#' split variable is categorical then this component contains the index of
#' \code{gbm_fit_obj$c.split} that describes the categorical split. If the node is a
#' terminal node then this is the prediction.} \item{LeftNode}{the index of the
#' row corresponding to the left node.} \item{RightNode}{the index of the row
#' corresponding to the right node.} \item{ErrorReduction}{the reduction in the
#' loss function as a result of splitting this node.} \item{Weight}{the total
#' weight of observations in the node. If weights are all equal to 1 then this
#' is the number of observations in the node.}
#' 
#' @author Greg Ridgeway \email{gregridgeway@@gmail.com}
#' @seealso \code{\link{gbm}}, \code{\link{gbm.object}}
#' @keywords print
#' @export 
#' 
pretty.GBMFit <- function(gbm_fit_obj, tree_index=1)
{
  # Initial checks
  check_if_natural_number(tree_index)
  if(tree_index >length(gbm_fit_obj$trees)) {
    stop("tree_index is out of range. Must be less than ", length(gbm_fit_obj$trees))
  }
  
  # Convert selected tree to data.frame
  tree <- data.frame(gbm_fit_obj$trees[[tree_index]])
  names(tree) <- c("SplitVar","SplitCodePred","LeftNode",
                   "RightNode","MissingNode","ErrorReduction",
                   "Weight","Prediction")
  row.names(tree) <- 0:(nrow(tree)-1)
  return(tree)
}