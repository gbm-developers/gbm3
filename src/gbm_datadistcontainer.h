//------------------------------------------------------------------------------
//
//  File:       gbmDataContainer.h
//
//  Description:   Header file for a class that stores the data and dist for
//  gbm.
//
//------------------------------------------------------------------------------

#ifndef GBMDATACONTAINER_H
#define GBMDATACONTAINER_H

//------------------------------
// Includes
//------------------------------
#include "datadistparams.h"
#include "dataset.h"
#include "databag.h"
#include "distribution.h"
#include "distribution_factory.h"
#include "tree.h"
#include "treeparams.h"
#include <Rcpp.h>
#include <vector>
#include <memory>

//------------------------------
// Class definition
//------------------------------
class CGBMDataDistContainer {
 public:
  //----------------------
  // Public Constructors
  //----------------------
  CGBMDataDistContainer(DataDistParams& datadist_config);

  //---------------------
  // Public destructor
  //---------------------
  ~CGBMDataDistContainer() {};

  //---------------------
  // Public Functions
  //---------------------
  void Initialize();
  double InitialFunctionEstimate();
  void ComputeResiduals(const double* kFuncEstimate, std::vector<double>& residuals);
  void ComputeBestTermNodePreds(const double* kFuncEstimate, std::vector<double>& residuals,
                                CCARTTree& tree);
  double ComputeDeviance(const double* kFuncEstimate,
                         bool is_validationset = false);
  double ComputeBagImprovement(const double* kFuncEstimate,
                               const double kShrinkage,
                               const std::vector<double>& kDeltaEstimate);
  void BagData();

  std::auto_ptr<CDistribution>& get_dist() { return distptr_; }
  const CDataset& get_data() const { return data_; }
  CDataset& get_data() { return data_; }
  const Bag& get_bag() const{ return databag_; }
  Bag& get_bag() { return databag_; }

 private:
  //-------------------
  // Private functions
  //-------------------
  void shift_datadist_to_train() {
    data_.shift_to_train();
    distptr_->ShiftDistPtrs(-(data_.get_trainsize()));
  }
  void shift_datadist_to_validation() {
    data_.shift_to_validation();
    distptr_->ShiftDistPtrs(data_.get_trainsize());
  }

  //-------------------
  // Private Variables
  //-------------------
  CDataset data_;
  Bag databag_;
  std::auto_ptr<DistributionFactory> distfactory_;
  std::auto_ptr<CDistribution> distptr_;
};

#endif  // GBMDATACONTAINER_H
