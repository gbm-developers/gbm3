//------------------------------------------------------------------------------
//  File:       quantile.h
//
//  Contents:   quantile regression for GBM.
//
//  History:    10/8/2006   Created by Brian Kriegler (bk@stat.ucla.edu)
//              6/11/2007   gregr merged with official gbm
//
//------------------------------------------------------------------------------

#ifndef QUANTILE_H
#define QUANTILE_H

//------------------------------
// Includes
//------------------------------
#include "distribution.h"
#include "locationm.h"
#include <algorithm>
#include <memory>

//------------------------------
// Class definition
//------------------------------
class CQuantile : public CDistribution {
 public:
  //---------------------
  // Factory Function
  //---------------------
  static CDistribution* Create(DataDistParams& distparams);

  //---------------------
  // Public destructor
  //---------------------
  virtual ~CQuantile();

  //---------------------
  // Public Functions
  //---------------------
  void ComputeWorkingResponse(const CDataset& kData, const Bag& kBag,
                              const double* kFuncEstimate,
                              std::vector<double>& residuals);

  double InitF(const CDataset& kData);

  void FitBestConstant(const CDataset& kData, const Bag& kBag,
                       const double* kFuncEstimate,
                       unsigned long num_terminalnodes,
                       std::vector<double>& residuals, CCARTTree& tree);

  double Deviance(const CDataset& kData, const Bag& kBag,
                  const double* kFuncEstimate);

  double BagImprovement(const CDataset& kData, const Bag& kBag,
                        const double* kFuncEstimate, const double kShrinkage,
                        const std::vector<double>& kDeltaEstimate);

 private:
  //----------------------
  // Private Constructors
  //----------------------
  CQuantile(double alpha);

  //-------------------
  // Private Variables
  //-------------------
  vector<double> vecd_;
  double alpha_;
  CLocationM mplocm_;
};

#endif  // QUANTILE_H
