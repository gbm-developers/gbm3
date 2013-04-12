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
gbmCluster(n)
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
  \item{id}{The interaction depth}
  \item{w}{The weights}
  \item{n}{The number of cores to use in the cluster.}
  \item{o}{The offset}
%  \item{verbose}{Whether or not to print output to screen}
%  \item{X, var.monotone, n.trees, n.minobsinnode, shrinkage, bag.fraction,
%        var.names, response.name, cv.group}{Arguments passed
%        through to gbm.fit}
} % Close \arguments
\details{
  These are functions used internally by \code{gbm} and not intended for
  direct use by the user.
}
