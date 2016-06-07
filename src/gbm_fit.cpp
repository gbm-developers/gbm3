//-----------------------------------
//
// File: gbm_fit
//
// Description: Object storing fitted Gbm model.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "gbm_fit.h"


//-----------------------------------
// Public Functions
//-----------------------------------
//-----------------------------------
// Function: GbmFit
//
// Returns: None
//
// Description: Constructor
//
// Parameters:
//-----------------------------------
GbmFit::GbmFit(const int kNumDataRows, const double kInitEstimate,
		  const int kNumTrees,
		  const Rcpp::NumericVector kPrevFuncEstimate):
		  func_estimate_(kNumDataRows),
		  training_errors_(kNumTrees, 0.0),
		  validation_errors_(kNumTrees, 0.0),
		  outofbag_improvement_(kNumTrees, 0.0),
		  set_of_trees_(kNumTrees),
		  initial_estimate_(kInitEstimate)
		  {
	if (ISNA(kPrevFuncEstimate[0]))  // check for old predictions
	{
		// set the initial value of F as a constant
		func_estimate_.fill(initial_estimate_);
	} else {
	if (kPrevFuncEstimate.size() != func_estimate_.size()) {
	  throw gbm_exception::InvalidArgument(
		  "old predictions are the wrong shape");
	}
	std::copy(kPrevFuncEstimate.begin(), kPrevFuncEstimate.end(), func_estimate_.begin());
	}
}


void GbmFit::AccumulateErrors(const int kTree, CGBMEngine& gbm) {
	std::vector<double> metrics = gbm.FitLearner(func_estimate_.begin());
	training_errors_[kTree] += metrics[0];
	validation_errors_[kTree] += metrics[1];
	outofbag_improvement_[kTree] += metrics[2];
}

void GbmFit::CreateTrees(const int kTree,  const int kCatSplitsOld, const CGBMEngine& gbm) {

	// Vectors defining a tree
	Rcpp::IntegerVector split_vars(gbm.size_of_fitted_tree());
	Rcpp::NumericVector split_values(gbm.size_of_fitted_tree());
	Rcpp::IntegerVector left_nodes(gbm.size_of_fitted_tree());
	Rcpp::IntegerVector right_nodes(gbm.size_of_fitted_tree());
	Rcpp::IntegerVector missing_nodes(gbm.size_of_fitted_tree());
	Rcpp::NumericVector error_reduction(gbm.size_of_fitted_tree());
	Rcpp::NumericVector weights(gbm.size_of_fitted_tree());
	Rcpp::NumericVector node_predictions(gbm.size_of_fitted_tree());

	gbm.GbmTransferTreeToRList(
	        split_vars.begin(), split_values.begin(), left_nodes.begin(),
	        right_nodes.begin(), missing_nodes.begin(), error_reduction.begin(),
	        weights.begin(), node_predictions.begin(), split_codes_, kCatSplitsOld);

	set_of_trees_[kTree] = Rcpp::List::create(
	        split_vars, split_values, left_nodes, right_nodes, missing_nodes,
	        error_reduction, weights, node_predictions);
}
