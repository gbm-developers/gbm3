//------------------------------------------------------------------------------
//
//  File:       poisson.h
//
//  Description:   poisson object for GBM.
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef POISSON_H
#define POISSON_H

//------------------------------
// Includes
//------------------------------
#include "distribution.h"
#include <Rmath.h>
#include <memory>

//------------------------------
// Class definition
//------------------------------
class CPoisson : public CDistribution {
 public:
  //---------------------
  // Factory Function
  //---------------------
  static CDistribution* Create(DataDistParams& distparams);

  //---------------------
  // Public destructor
  //---------------------
  virtual ~CPoisson();

  //---------------------
  // Public Functions
  //---------------------
  void ComputeWorkingResponse(const CDataset& kData,
		  	  	  	  	  	  const Bag& kBag,
                              const double* kFuncEstimate, std::vector<double>& residual);

  double Deviance(const CDataset& kData, const Bag& kBag, const double* kFuncEstimate);

  double InitF(const CDataset& kData);

  void FitBestConstant(const CDataset& kData, const Bag& kBag, const double* kFuncEstimate,
                       unsigned long num_terminalnodes, std::vector<double>& residual,
                       CCARTTree& tree);

  double BagImprovement(const CDataset& kData, const Bag& kBag, const double* kFuncEstimate,
                        const double kShrinkage,
                        const std::vector<double>& kDeltaEstimate);

 private:
  //----------------------
  // Private Constructors
  //----------------------
  CPoisson();
};

#endif  // POISSON_H
