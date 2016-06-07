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
  void AccumulateErrors(const int kTree, CGBMEngine& gbm);
  void CreateTrees(const int kTree, const int kCatSplitsOld, const CGBMEngine& gbm);
  Rcpp::List ROutput() {
	return  Rcpp::List::create(
	       _["initF"] = initial_estimate_, _["fit"] = func_estimate_,
	       _["train.error"] = training_errors_, _["valid.error"] = validation_errors_,
	       _["oobag.improve"] = outofbag_improvement_, _["trees"] = set_of_trees_,
	       _["c.splits"] = split_codes_);
  }

  // Inlined Getters
  double get_tree_training_error(const int kTree) const {
	  return training_errors_[kTree];
  }
  double get_tree_valid_error(const int kTree) const {
	  return validation_errors_[kTree];
  }
  double get_tree_oobag_improv(const int kTree) const {
	  return outofbag_improvement_[kTree];
  }


 private:
  //---------------------
  // Private Variables
  //---------------------
  VecOfVectorCategories split_codes_;
  Rcpp::NumericVector func_estimate_; // Fitted function
  Rcpp::NumericVector training_errors_;
  Rcpp::NumericVector validation_errors_;
  Rcpp::NumericVector outofbag_improvement_;
  Rcpp::GenericVector set_of_trees_;
  double initial_estimate_;
};

#endif  // GBMFIT_H
