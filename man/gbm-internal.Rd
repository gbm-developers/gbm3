\name{gbm-internal}
\alias{guessDist}
\alias{getStratify}
\alias{getCVgroup}
\alias{checkMissing}
\alias{checkID}
\alias{checkWeights}
\alias{checkOffset}
\alias{getVarNames}
\alias{gbmCluster}
\alias{gbmDoFold}
\title{gbm internal functions}
\description{Helper functions for preprocessing data prior to
             building the model}
\usage{
guessDist(y)
getCVgroup(distribution, class.stratify.cv, y, i.train, cv.folds, group)
getStratify(strat, d)
checkMissing(x, y)
checkWeights(w, n)
checkID(id)
checkOffset(o, y)
getVarNames(x)
gbmCluster(n, ncv, verbose)
gbmDoFold(X, i.train, x, y, offset, distribution, w, var.monotone, n.trees,
         interaction.depth, n.minobsinnode, shrinkage, bag.fraction,
         cv.group, var.names, response.name, group)
}
\arguments{
  \item{y}{The response variable}
  \item{d, distribution}{The distribution, either specified by the user or
        implied}
  \item{class.stratify.cv}{Whether or not to stratify, if provided by
        the user}
  \item{i.train}{Computed internally by \code{gbm}}
  \item{group}{The group, if using \code{distibution='pairwise'}}
  \item{strat}{Whether or not to stratify}
  \item{cv.folds}{The number of cross-validation folds}
  \item{x}{The design matrix}
  \item{id, interaction.depth}{The interaction depth}
  \item{w}{The weights}
  \item{o, offset}{The offset}
  \item{n, ncv}{The length of the response variable, or the number of
    cores and the number of cross-validations folds.}
  \item{verbose}{Whether or not to print output to screen}
  \item{X, var.monotone, n.trees, n.minobsinnode, shrinkage, bag.fraction,
        var.names, response.name, cv.group}{Arguments passed
        through to gbm.fit}
} % Close \arguments
\details{
  These are functions used internally by \code{gbm} and not intended for
  direct use by the user.
}
