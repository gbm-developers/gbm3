//------------------------------------------------------------------------------
//
//  File:       gbm.h
//
//  License:    GNU GPL (version 2 or later)
//
//  Description:   Entry point for gbm.dll
//
//  Owner:      gregr@rand.org
//
//  History:    2/14/2003   gregr created
//              6/11/2007   gregr added quantile regression
//                          written by Brian Kriegler
//
//------------------------------------------------------------------------------

#ifndef __gbm_h__
#define __gbm_h__

//------------------------------
// Includes
//------------------------------
#include <vector>
#include <string>
#include <memory>
#include <Rcpp.h>
#include "dataset.h"
#include "distribution.h"
#include "gbmTreeComps.h"
#include "bernoulli.h"
#include "adaboost.h"
#include "poisson.h"
#include "gaussian.h"
#include "coxph.h"
#include "laplace.h"
#include "quantile.h"
#include "tdist.h"
#include "pairwise.h"
#include "gbm_engine.h"
#include "locationm.h"
#include "huberized.h"
#include "gamma.h"
#include "tweedie.h"

CDistribution* gbm_setup
(
    const CDataset& data,
    SEXP radMisc,
    const std::string& family,
    int cTrain,
    int& cGroups
);


#endif // __gbm_h__
