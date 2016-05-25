//------------------------------------------------------------------------------
//
//  File:       huberized.h
//
//  Description:   huberized hinge loss object for GBM.
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//------------------------------------------------------------------------------

#ifndef HUBERIZED_H
#define HUBERIZED_H

//------------------------------
// Includes
//------------------------------
#include "distribution.h"
#include <memory>

//------------------------------
// Class definition
//------------------------------
class CHuberized : public CDistribution {
 public:
  //---------------------
  // Factory Function
  //---------------------
  static CDistribution* Create(DataDistParams& distparams);

  //---------------------
  // Public destructor
  //---------------------
  virtual ~CHuberized();

  //---------------------
  // Public Functions
  //---------------------
  void ComputeWorkingResponse(const CDataset& kData,
                              const double* kFuncEstimate, double* residuals);

  double Deviance(const CDataset& kData, const double* kFuncEstimate);

  double InitF(const CDataset& kData);

  void FitBestConstant(const CDataset& kData, const double* kFuncEstimate,
                       unsigned long num_terminalnodes, double* residuals,
                       CCARTTree& tree);

  double BagImprovement(const CDataset& kData, const double* kFuncEstimates,
                        const double kShrinkage, const double* kDeltaEstimates);

 private:
  //----------------------
  // Private Constructors
  //----------------------
  CHuberized();
};

#endif  // HUBERIZED_H
