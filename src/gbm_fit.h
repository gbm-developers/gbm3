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

using Rcpp::_;

//------------------------------
// Class definition
//------------------------------
class GbmFit {
 public:
  //----------------------
  // Public Constructors
  //----------------------
  GbmFit(const int kNumDataRows, const double kInitEstimate,
         const int kNumTrees, const Rcpp::NumericVector& kPrevFuncEstimate);

  //---------------------
  // Public Functions
  //---------------------
  void accumulate(CGBMEngine& gbm);
  void CreateTreeRepresentation(const int kCatSplitsOld);
  Rcpp::List ROutput();

  // Inlined functions
  double get_tree_training_error() { return training_errors_[tree_count_]; }
  double get_tree_valid_error() { return validation_errors_[tree_count_]; }
  double get_tree_oobag_improv() { return outofbag_improvement_[tree_count_]; }
  void increment_count() { tree_count_++; }

 private:
  //---------------------
  // Private Variables
  //---------------------
  VecOfVectorCategories split_codes_;
  std::auto_ptr<FittedLearner> current_fit_;
  Rcpp::NumericVector training_errors_;
  Rcpp::NumericVector validation_errors_;
  Rcpp::NumericVector outofbag_improvement_;
  Rcpp::NumericVector func_estimate_;  // Fitted function
  Rcpp::GenericVector set_of_trees_;
  double initial_estimate_;
  unsigned long tree_count_;
};

#endif  // GBMFIT_H
