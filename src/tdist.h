//------------------------------------------------------------------------------
//
//  File:       tdist.h
//
//  Description:   Distribution object to implement t-distribution
//
//  History:    04/04/2008   Created
//
//------------------------------------------------------------------------------

#ifndef TDIST_H
#define TDIST_H

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
class CTDist : public CDistribution {
 public:
  //---------------------
  // Factory Function
  //---------------------
  static CDistribution* Create(DataDistParams& distarams);

  //---------------------
  // Public destructor
  //---------------------
  virtual ~CTDist();

  //---------------------
  // Public Functions
  //---------------------
  void ComputeWorkingResponse(const CDataset& kData,
                              const double* kFuncEstimate, std::vector<double>& residuals);

  double InitF(const CDataset& kData);

  void FitBestConstant(const CDataset& kData, const double* kFuncEstimate,
                       unsigned long num_terminalnodes, std::vector<double>& residuals,
                       CCARTTree& tree);

  double Deviance(const CDataset& kData, const double* kFuncEstimates);

  double BagImprovement(const CDataset& kData, const double* kFuncEstimate,
                        const double kShrinkage, const std::vector<double>& kDeltaEstimate);

 private:
  //----------------------
  // Private Constructors
  //----------------------
  CTDist(double nu);

  //-------------------
  // Private Variables
  //-------------------
  double m_nu_;
  CLocationM mplocm_;
};

#endif  // TDIST_H
