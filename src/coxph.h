//------------------------------------------------------------------------------
//
//  File:       coxph.h
//
//  Description:   Cox proportional hazard object
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef COXPH_H
#define COXPH_H

//------------------------------
// Includes
//------------------------------

#include "distribution.h"
#include <memory>

//------------------------------
// Class Forwards and Enums
//------------------------------
class GenericCoxState;

//------------------------------
// Class definition
//------------------------------
class CCoxPH : public CDistribution {
 public:
  //---------------------
  // Factory Function
  //---------------------
  static CDistribution* Create(DataDistParams& distparams);

  //---------------------
  // Public destructor
  //---------------------
  ~CCoxPH(){};

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

  int TieApproxMethod() const { return tiedtimesmethod_; }
  double PriorCoeffVar() const { return kPriorCoeffVariation_; }

 private:
  //----------------------
  // Private Constructors
  //----------------------
  CCoxPH(bool is_startstop, int tiesmethod, double priorcoeff);


  //-------------------
  // Private Variables
  //-------------------
  const bool kStartStopCase_;
  const double kPriorCoeffVariation_;
  int tiedtimesmethod_;
  std::auto_ptr<GenericCoxState> coxstate_methods_;
};

#endif  // COXPH_H
