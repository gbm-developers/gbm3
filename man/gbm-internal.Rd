\name{gbm-internal}
\alias{guessDist}
\alias{getStratify}
\alias{getCVgroup}
\title{gbm internal functions}
\description{Helper functions for preprocessing data prior to
             building the model}
\usage{
guessDist(y)
getCVgroup(distribution, class.stratify.cv, y, i.train, cv.folds, group)
getStratify(strat, d)
}
\arguments{
  \item{y}{The response variable}
  \item{distribution}{The distribution, either specified by the user or
    implied}
  \item{class.stratify.cv}{Whether or not to stratify, if provided by
    the user}
  \item{i.train}{Computed internally by \code{gbm}}
  \item{group}{The group, if using \code{distibution='pairwise'}}
  \item{d}{The distribution}
  \item{strat}{Whether or not to stratify}
  \item{cv.folds}{The number of cross-validation folds}
} % Close \arguments
\details{
  These are functions used internally by \code{gbm} and not intended for
  direct use by the user.
}
