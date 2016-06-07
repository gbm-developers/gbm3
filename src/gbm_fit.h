//------------------------------------------------------------------------------
//
//  File:       gbm_fit.h
//
//  Description:   header for gbm fit output class.
//
//  Owner: James Hickey
//
//  History:    07/06/2016  jhickey: created
//
//------------------------------------------------------------------------------

#ifndef GBMFIT_H
#define GBMFIT_H

//------------------------------
// Includes
//------------------------------
#include "gbm_engine.h"
#include <Rcpp.h>

using Rcpp::_ ;

//------------------------------
// Class definition
//------------------------------
class GbmFit {
 public:
  //----------------------
  // Public Constructors
  //----------------------
  GbmFit(const int kNumDataRows, const double kInitEstimate,
		  const int kNumTrees,
		  const Rcpp::NumericVector kPrevFuncEstimate);

  //---------------------
  // Public destructor
  //---------------------
  ~GbmFit() {};

  //---------------------
  // Public Functions
  //---------------------
  void accumulate(CGBMEngine& gbm);
  void CreateTreeRepresentation(const int kCatSplitsOld);
  Rcpp::List ROutput();

  // Inlined functions
  double get_tree_training_error() const {
	  return current_fit_[tree_count_].training_error;
  }
  double get_tree_valid_error() const {
	  return current_fit_[tree_count_].validation_error;
  }
  double get_tree_oobag_improv() const {
	  return current_fit_[tree_count_].oobag_improvement;
  }
  void increment_count() { tree_count_++; }


 private:
  //---------------------
  // Private Variables
  //---------------------
  VecOfVectorCategories split_codes_;
  std::vector<FitStruct> current_fit_;
  Rcpp::NumericVector func_estimate_; // Fitted function
  Rcpp::GenericVector set_of_trees_;
  double initial_estimate_;
  unsigned long tree_count_;
};

#endif  // GBMFIT_H
