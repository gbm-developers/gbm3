//------------------------------------------------------------------------------
//
//  File:       genericCoxState.h
//
//  Description: abstract class defining the generic methods associated with the
//  CoxPh model.
//
//	Author: 	James Hickey
//------------------------------------------------------------------------------

#ifndef GENERICCOXSTATE_H
#define GENERICCOXSTATE_H

//-----------------------------------
// Definitions
//-----------------------------------
#define frac .00000001
#define recenter 50

//------------------------------
// Includes
//------------------------------
#include "dataset.h"
#include <Rcpp.h>

//------------------------------
// Generic Dispatch Definition
//------------------------------
class GenericCoxState {
 public:
  //----------------------
  // Public Constructors
  //----------------------
  GenericCoxState(){};

  //---------------------
  // Public destructor
  //---------------------
  virtual ~GenericCoxState(){};

  //---------------------
  // Public Functions
  //---------------------
  virtual void ComputeWorkingResponse(const CDataset& kData,
                                      const double* kFuncEstimate,
                                      double* residuals) = 0;

  virtual void FitBestConstant(const CDataset& kData,
                               const double* kFuncEstimate,
                               unsigned long num_terminalnodes,
                               double* residuals, CCARTTree& tree) = 0;

  virtual double Deviance(const long kNumRowsInSet, const CDataset& kData,
                          const double* kFuncEstimate) = 0;

  virtual double BagImprovement(const CDataset& kData,
                                const double* kFuncEstimate,
                                const double kShrinkage,
                                const double* kDeltaEstimate) = 0;
};
#endif  // GENERICCOXSTATE_H
