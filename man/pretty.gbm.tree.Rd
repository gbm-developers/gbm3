\name{pretty.gbm.tree}
\alias{pretty.gbm.tree}
\title{ Print gbm tree components }
\description{
\code{gbm} stores the collection of trees used to construct the model in a 
compact matrix structure. This function extracts the information from a single
tree and displays it in a slightly more readable form. This function is mostly
for debugging purposes and to satisfy some users' curiosity.
}
\usage{
pretty.gbm.tree(object, i.tree = 1)
}
\arguments{
  \item{object}{ a \code{\link{gbm.object}} initially fit using \code{\link{gbm}}}
  \item{i.tree}{ the index of the tree component to extract from \code{object} 
                 and display }
}
\value{
\code{pretty.gbm.tree} returns a data frame. Each row corresponds to a node in 
the tree. Columns indicate
\item{SplitVar}{index of which variable is used to split. -1 indicates a 
terminal node.}
\item{SplitCodePred}{if the split variable is continuous then this component
is the split point. If the split variable is categorical then this component
contains the index of \code{object$c.split} that describes the categorical
split. If the node is a terminal node then this is the prediction.}
\item{LeftNode}{the index of the row corresponding to the left node.}
\item{RightNode}{the index of the row corresponding to the right node.}
\item{ErrorReduction}{the reduction in the loss function as a result of 
splitting this node.}
\item{Weight}{the total weight of observations in the node. If weights are all
equal to 1 then this is the number of observations in the node.}
}
\author{Greg Ridgeway \email{gregridgeway@gmail.com}}
\seealso{
\code{\link{gbm}}, \code{\link{gbm.object}}
}
\keyword{ print }
