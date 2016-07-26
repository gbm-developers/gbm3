//------------------------------------------------------------------------------
//
//  File:       gaussian.h
//
//  Contents:   gaussian distribution for GBM.
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//------------------------------------------------------------------------------

#ifndef GAUSSIAN_H
#define GAUSSIAN_H

//------------------------------
// Includes
//------------------------------
#include "distribution.h"
#include <memory>

//------------------------------
// Class definition
//------------------------------
class CGaussian : public CDistribution {
 public:
  //---------------------
  // Factory Function
  //---------------------
  static CDistribution* Create(DataDistParams& distparams);

  //---------------------
  // Public destructor
  //---------------------
  virtual ~CGaussian();

  //---------------------
  // Public Functions
  //---------------------
  void ComputeWorkingResponse(const CDataset& kData, const Bag& kBag,
                              const double* kFuncEstimate,
                              std::vector<double>& residuals);

  double InitF(const CDataset& data);

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
  CGaussian(const parallel_details& parallel);
};

#endif  // GAUSSIAN_H
