//------------------------------------------------------------------------------
//
//  File:       gamma.h
//
//  Description: gamma distribution
//
//------------------------------------------------------------------------------

#ifndef GAMMA_H
#define GAMMA_H

//------------------------------
// Includes
//------------------------------
#include "distribution.h"
#include <Rmath.h>
#include <memory>

//------------------------------
// Class definition
//------------------------------
class CGamma : public CDistribution {
 public:
  //---------------------
  // Factory Function
  //---------------------
  static CDistribution* Create(DataDistParams& distparams);

  //---------------------
  // Public destructor
  //---------------------
  virtual ~CGamma();

  //---------------------
  // Public Functions
  //---------------------
  void ComputeWorkingResponse(const CDataset& kData,
                              const double* kFuncEstimate, std::vector<double>& residuals);

  double InitF(const CDataset& kData);

  void FitBestConstant(const CDataset& Data, const double* kFuncEstimate,
                       unsigned long num_terminalnodes, std::vector<double>& residuals,
                       CCARTTree& tree);

  double Deviance(const CDataset& kData, const double* kFuncEstimate);

  double BagImprovement(const CDataset& kData, const double* kFuncEstimate,
                        const double kShrinkage, const std::vector<double>& kDeltaEstimate);

 private:
  //----------------------
  // Private Constructors
  //----------------------
  CGamma();
};

#endif  // GAMMA_H
