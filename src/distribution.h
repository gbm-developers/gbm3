//------------------------------------------------------------------------------
//
//  File:       distribution.h
//
//  Description: Header file for distribution object used in GBM. This class
//    is the parent class from which all other distributions inherit.
//
//
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//				22/02/2016  jhickey: modified to implement
// factory
// pattern
//------------------------------------------------------------------------------

#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

//------------------------------
// Includes
//------------------------------
#include "dataset.h"
#include "databag.h"
#include "node.h"
#include "parallel_details.h"
#include "tree.h"
#include <vector>
#include <Rcpp.h>

//------------------------------
// Class definition
//------------------------------
class CDistribution {
 public:
  //----------------------
  // Public Constructors
  //----------------------
  CDistribution();
  CDistribution(const parallel_details& parallel);

  //---------------------
  // Public destructor
  //---------------------
  virtual ~CDistribution();

  //---------------------
  // Public Functions
  //---------------------
  int get_num_groups() const { return num_groups_; };
  void set_num_groups(int groupval) { num_groups_ = groupval; };

  // shifts the ptrs() as appropriate
  template <typename T>
  T* shift_ptr(T* x, std::ptrdiff_t y) {
    if (x) {
      return x + y;
    } else {
      return x;
    }
  }

  //---------------------
  // Public Virtual Functions
  //---------------------
  virtual void Initialize(const CDataset& kData) {

    // Set up multi map
    for (unsigned long i = 0;
         i < (kData.get_trainsize()); i++) {
      obsid_to_row_.insert(pair<int, int>(kData.get_row_observation_id(i), i));
    }
  };
  virtual void ComputeWorkingResponse(const CDataset& kData, const Bag& kBag,
                                      const double* kFuncEstimate,
                                      std::vector<double>& residuals) = 0;

  virtual double InitF(const CDataset& kData) = 0;

  virtual double Deviance(const CDataset& kData, const Bag& kBag,
                          const double* kFuncEstimate) = 0;

  virtual void FitBestConstant(const CDataset& kData, const Bag& kBag,
                               const double* kFuncEstimate,
                               unsigned long num_terminalnodes,
                               std::vector<double>& residuals,
                               CCARTTree& tree) = 0;

  virtual double BagImprovement(
      const CDataset& kData, const Bag& kBag, const double* kFuncEstimate,
      const double kShrinkage,
      const std::vector<double>& kDeltaFuncEstimate) = 0;

  virtual void BagData(const CDataset& kData, Bag& bag);
  virtual void ShiftDistPtrs(unsigned long shift){};

  int get_num_threads() const { return parallel_.get_num_threads(); }
  int get_array_chunk_size() const { return parallel_.get_array_chunk_size(); }
  
 private:
  //---------------------
  // Private Variables
  //---------------------
  parallel_details parallel_;
  int num_groups_;
  std::multimap<int, int> obsid_to_row_;  // Map from observation unit to row
};

#endif  // DISTRIBUTION_H
