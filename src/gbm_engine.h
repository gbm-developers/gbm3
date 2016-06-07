//------------------------------------------------------------------------------
//
//  File:       gbm_engine.h
//
//  Description:   Header file for Gradient Boosting Engine.
//
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef GBMENGINE_H
#define GBMENGINE_H

//------------------------------
// Includes
//------------------------------
#include "datadistparams.h"
#include "fitstruct.h"
#include "gbm_datadistcontainer.h"
#include "tree.h"
#include "treeparams.h"
#include <memory>
#include <Rcpp.h>
#include <vector>

//------------------------------
// Class definition
//------------------------------
class CGBMEngine {
 public:
  //----------------------
  // Public Constructors
  //----------------------
  CGBMEngine(DataDistParams& datadistparams, TreeParams& treeparams);

  //---------------------
  // Public destructor
  //---------------------
  ~CGBMEngine();

  //---------------------
  // Public Functions
  //---------------------
  FitStruct FitLearner(double* func_estimate);
  double initial_function_estimate() {
    return datacontainer_.InitialFunctionEstimate();
  };

 private:
  //-------------------
  // Private Variables
  //-------------------
  CGBMDataDistContainer datacontainer_;
  TreeParams& tree_params_;

  // Residuals and adjustments to function estimate
  std::vector<double> residuals_;
};

#endif  // GBMENGINE_H
