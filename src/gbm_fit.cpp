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
		  current_fit_(kNumTrees, FitStruct()),
		  func_estimate_(kNumDataRows),
		  set_of_trees_(kNumTrees),
		  initial_estimate_(kInitEstimate),
		  tree_count_(0)
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


void GbmFit::accumulate(CGBMEngine& gbm) {
	current_fit_[tree_count_] = gbm.FitLearner(func_estimate_.begin());
}

void GbmFit::CreateTreeRepresentation(const int kCatSplitsOld) {

    // Vectors defining a tree
	Rcpp::IntegerVector split_vars(current_fit_[tree_count_].fitted_tree->size_of_tree());
	Rcpp::NumericVector split_values(current_fit_[tree_count_].fitted_tree->size_of_tree());
	Rcpp::IntegerVector left_nodes(current_fit_[tree_count_].fitted_tree->size_of_tree());
	Rcpp::IntegerVector right_nodes(current_fit_[tree_count_].fitted_tree->size_of_tree());
	Rcpp::IntegerVector missing_nodes(current_fit_[tree_count_].fitted_tree->size_of_tree());
	Rcpp::NumericVector error_reduction(current_fit_[tree_count_].fitted_tree->size_of_tree());
	Rcpp::NumericVector weights(current_fit_[tree_count_].fitted_tree->size_of_tree());
	Rcpp::NumericVector node_predictions(current_fit_[tree_count_].fitted_tree->size_of_tree());

	current_fit_[tree_count_].fitted_tree->TransferTreeToRList(*current_fit_[tree_count_].data_for_fit,
			split_vars.begin(), split_values.begin(), left_nodes.begin(),
			right_nodes.begin(), missing_nodes.begin(), error_reduction.begin(),
			weights.begin(), node_predictions.begin(), split_codes_, kCatSplitsOld);

	set_of_trees_[tree_count_] = Rcpp::List::create(
			split_vars, split_values, left_nodes, right_nodes, missing_nodes,
			error_reduction, weights, node_predictions);
}

Rcpp::List GbmFit::ROutput() {
	Rcpp::NumericVector training_errors(tree_count_, 0.0);
    Rcpp::NumericVector validation_errors(tree_count_, 0.0);
	Rcpp::NumericVector outofbag_improvement(tree_count_, 0.0);

	for(unsigned long treenum = 0; treenum < tree_count_; treenum++) {
		training_errors[treenum] += current_fit_[treenum].training_error;
		validation_errors[treenum] += current_fit_[treenum].validation_error;
		outofbag_improvement[treenum] += current_fit_[treenum].oobag_improvement;
	}

	return  Rcpp::List::create(
	       _["initF"] = initial_estimate_, _["fit"] = func_estimate_,
	       _["train.error"] = training_errors, _["valid.error"] = validation_errors,
	       _["oobag.improve"] = outofbag_improvement, _["trees"] = set_of_trees_,
	       _["c.splits"] = split_codes_);
 }
