// GBM by Greg Ridgeway  Copyright (C) 2003

#include "datadistparams.h"
#include "gbm_engine.h"
#include "gbm_fit.h"
#include "parallel_details.h"
#include "treeparams.h"
#include <memory>
#include <utility>
#include <Rcpp.h>

namespace {
class NodeStack {
 public:
  bool empty() const { return stack.empty(); }

  std::pair<int, double> pop() {
    std::pair<int, double> result = stack.back();
    stack.pop_back();
    return result;
  }

  void push(int nodeIndex, double weight = 0.0) {
    stack.push_back(std::make_pair(nodeIndex, weight));
  }

 private:
  std::vector<std::pair<int, double> > stack;
};


  inline parallel_details parallel_details_wrap(SEXP src) {
    Rcpp::List details(src);
    int num_threads = details["n.threads"];
    int array_chunk_size = details["arrayChunkSize"];
    return parallel_details(num_threads, array_chunk_size);
  }

}

//----------------------------------------
// Functions- Public
//----------------------------------------
extern "C" {

//-----------------------------------
// Function: gbm
//
// Returns: R List containing: the initial function estimate, the fit,
//			training errors, validation errors, out of bag
//			improvement, the trees and the categorical
//			splits.
//
// Description: Fits a gbm model to data.
//
// Parameters:
//  response  - SEXP containing the response of each data-point - accessed via
//				double ptr
//  offset_vec - SEXP containing the offset applied to each response - accessed
//  via
//				double ptr - NA for no offset
//  covariates - SEXP containing the predictor values - becomes
//  Rcpp::NumericMatrix
//  covar_order - SEXP containing the order of predictor values to
//  			be used in GBM formula - accessed via int ptr.
//  sorted_vec - SEXP indicating the ordering of observations for CoxPH
//				 , this is stored an int ptr in CoxPH.
//  strata_vec - SEXP indicating which strata observations are in -
//				 this is used with CoxPH as an int ptr.
//  obs_weight  - SEXP containing weights to be used in fitting
//  				process - accessed via double ptr.
//  misc - SEXP list object containing distribution dependent data
//		   , depending on the distribution this may be a double or int.
//		   It is NA if there is no misc.
// prior_coeff_var - SEXP containing a double specifying a prior node
//					prediction value for the CoxPH model.
//  row_to_obs_id - SEXP storing integer ids mapping each row to
//          		its observation - becomes Rcpp::IntegerVector.
//  var_classes  - R object (SEXP) containing the variable classes,
//				  this is stored as a Rcpp::IntegerVector in
// dataset.cpp.
//  monotonicity_vec - SEXP containing +-1/0 indicating whether the
//  					variables are:
//  monotone increasing (+1), decreasing (-1) or arbitrary (0) with the response
// 					 variable. - accessed via int ptr
//  dist_family - SEXP specifying which distribution is used to model data -
//  string.
//  num_trees - SEXP specifying how many trees to fit in model - ulong.
//  tree_depth - SEXP defining maximum depth of each tree - ulong.
//  min_num_node_obs - SEXP specifying minimum num of obs a node must have -
//  ulong.
//  shrinkageconstant - SEXP defining the shrinkage applied to each tree fit -
//  double.
//  fraction_inbag - SEXP, the fraction of observations to be used in tree
//  fitting - double.
//	num_rows_in_training - SEXP containing ulong specify the number of data
// points in
//							training set
//  num_obs_in_training - SEXP containing ulong specifying number of unique
//  observations
//						  that make up training set
//  number_offeatures - SEXP containing ulong specifying number of features for
//  use in
//						tree growing.
//  prev_func_estimate - SEXP containing previous function estimate -const
//  Rcpp::NumericVector
//  prev_category_splits - SEXP containing const int specifying the last
//  categorical
//						 split code of a previous fit
// for
// further
// fitting - if
//						 first fit it is 0.
//  prev_trees_fitted -  SEXP containing const int - number of previous trees
//  fitted.
//  par_details - SEXP giving details about parallelization
//  isverbose - SEXP which is a const bool specifying whether the fitting should
//  be
//				silent or not.
//-----------------------------------
SEXP gbm(SEXP response, SEXP offset_vec, SEXP covariates, SEXP covar_order,
         SEXP sorted_vec, SEXP strata_vec, SEXP obs_weight, SEXP misc,
         SEXP prior_coeff_var, SEXP row_to_obs_id, SEXP var_classes,
         SEXP monotonicity_vec, SEXP dist_family, SEXP num_trees,
         SEXP tree_depth, SEXP min_num_node_obs, SEXP shrinkageconstant,
         SEXP fraction_inbag, SEXP num_rows_in_training,
         SEXP num_obs_in_training, SEXP number_offeatures,
         SEXP prev_func_estimate, SEXP prev_category_splits,
         SEXP prev_trees_fitted, SEXP par_details, SEXP isverbose) {
  BEGIN_RCPP

  // Set up consts for tree fitting and transfer to R API
  const int kNumTrees = Rcpp::as<int>(num_trees);
  const int kCatSplitsOld = Rcpp::as<int>(prev_category_splits);
  const int kTreesOld = Rcpp::as<int>(prev_trees_fitted);
  const bool kIsVerbose = Rcpp::as<bool>(isverbose);
  const Rcpp::NumericVector kPrevFuncEst(prev_func_estimate);

  // extract parallelization info in one place
  // as it's used by both the distribution and the tree
  const parallel_details parallel(parallel_details_wrap(par_details));

  // Set up parameters for initialization
  DataDistParams datadistparams(
      response, offset_vec, covariates, covar_order, sorted_vec, strata_vec,
      obs_weight, misc, prior_coeff_var, row_to_obs_id, var_classes,
      monotonicity_vec, dist_family, fraction_inbag, num_rows_in_training,
      num_obs_in_training, number_offeatures, parallel);
  TreeParams treeparams(tree_depth, min_num_node_obs, shrinkageconstant,
                        num_rows_in_training, parallel);
  Rcpp::RNGScope scope;

  // Initialize GBM engine
  CGBMEngine gbm(datadistparams, treeparams);

  // Initialize the output object
  GbmFit gbmfit(datadistparams.response.nrow(), gbm.initial_function_estimate(),
                kNumTrees, kPrevFuncEst);

  if (kIsVerbose) {
    Rprintf("Iter   TrainDeviance   ValidDeviance   StepSize   Improve\n");
  }
  for (int treenum = 0; treenum < kNumTrees; treenum++) {
    Rcpp::checkUserInterrupt();

    // Calculate Errors
    gbmfit.accumulate(gbm);

    // Create Trees
    gbmfit.CreateTreeRepresentation(kCatSplitsOld);

    // print the information
    if ((kIsVerbose) &&
        ((treenum <= 9) || (0 == (treenum + 1 + kTreesOld) % 20) ||
         (treenum == kNumTrees - 1))) {
      Rprintf("%6d %13.4f %15.4f %10.4f %9.4f\n", treenum + 1 + kTreesOld,
              gbmfit.get_tree_training_error(), gbmfit.get_tree_valid_error(),
              treeparams.shrinkage, gbmfit.get_tree_oobag_improv());
    }

    // Increment internal count
    gbmfit.increment_count();
  }
  if (kIsVerbose) Rprintf("\n");

  return gbmfit.ROutput();
  END_RCPP
}

//-----------------------------------
// Function: gbm_pred
//
// Returns: predicted function - a Rcpp::NumericVector
//
// Description: Makes predictions using a previously fit
//			gbm model and data.
//
// Parameters:
//  covariates - SEXP containing the predictor values - becomes
//				  const Rcpp::NumericMatrix.
//  num_trees - SEXP containing an int or vector of ints specifying the number
//  of
//				  trees to make predictions on - stored as const
// Rcpp::IntegerVector.
//  initial_func_est - SEXP specifying initial prediction estimate - double.
//  fitted_trees - SEXP containing lists defining the previously fitted trees -
//					stored as const Rcpp::GenericVector.
//  categorical_splits - SEXP containing  list of the categories of the split
//  variables
//						defining a tree - stored as a
// const
// Rcpp::GenericVector.
//  variable_type -  SEXP containing integers specifying whether the variable
//					is continuous/nominal- stored as const
// Rcpp::IntegerVector.
//  ret_single_tree_res - SEXP which is a const bool specifying whether the
//  fitting should be
//				only a single tree.
//-----------------------------------

SEXP gbm_pred(SEXP covariates, SEXP num_trees, SEXP initial_func_est,
              SEXP fitted_trees, SEXP categorical_splits, SEXP variable_type,
              SEXP ret_single_tree_res) {
  BEGIN_RCPP
  int tree_num = 0;
  int obs_num = 0;
  const Rcpp::NumericMatrix kCovarMat(covariates);
  const int kNumCovarRows = kCovarMat.nrow();
  const Rcpp::IntegerVector kTrees(num_trees);
  const Rcpp::GenericVector kFittedTrees(fitted_trees);
  const Rcpp::IntegerVector kVarType(variable_type);
  const Rcpp::GenericVector kSplits(categorical_splits);
  const bool kSingleTree = Rcpp::as<bool>(ret_single_tree_res);
  const int kPredIterations = kTrees.size();
  int prediction_iteration = 0;
  int class_num = 0;

  if ((kCovarMat.ncol() != kVarType.size())) {
    throw gbm_exception::InvalidArgument("shape mismatch");
  }

  Rcpp::NumericVector predicted_func(kNumCovarRows * kPredIterations);

  // initialize the predicted values
  if (!kSingleTree) {
    std::fill(predicted_func.begin(), predicted_func.begin() + kNumCovarRows,
              Rcpp::as<double>(initial_func_est));
  } else {
    predicted_func.fill(0.0);
  }
  tree_num = 0;
  for (prediction_iteration = 0; prediction_iteration < kTrees.size();
       prediction_iteration++) {
    const int kCurrTree = kTrees[prediction_iteration];
    if (kSingleTree) tree_num = kCurrTree - 1;
    if (!kSingleTree && (prediction_iteration > 0)) {
      // copy over from the last rcTrees
      std::copy(
          predicted_func.begin() + kNumCovarRows * (prediction_iteration - 1),
          predicted_func.begin() + kNumCovarRows * prediction_iteration,
          predicted_func.begin() + kNumCovarRows * prediction_iteration);
    }
    while (tree_num < kCurrTree) {
      const Rcpp::GenericVector kThisFitTree = kFittedTrees[tree_num];
      const Rcpp::IntegerVector kThisSplitVar = kThisFitTree[0];
      const Rcpp::NumericVector kThisSplitCode = kThisFitTree[1];
      const Rcpp::IntegerVector kThisLeftNode = kThisFitTree[2];
      const Rcpp::IntegerVector kThisRightNode = kThisFitTree[3];
      const Rcpp::IntegerVector kThisMissingNode = kThisFitTree[4];

      for (obs_num = 0; obs_num < kNumCovarRows; obs_num++) {
        int iCurrentNode = 0;
        while (kThisSplitVar[iCurrentNode] != -1) {
          const double dX =
              kCovarMat[kThisSplitVar[iCurrentNode] * kNumCovarRows + obs_num];
          // missing?
          if (ISNA(dX)) {
            iCurrentNode = kThisMissingNode[iCurrentNode];
          }
          // continuous?
          else if (kVarType[kThisSplitVar[iCurrentNode]] == 0) {
            if (dX < kThisSplitCode[iCurrentNode]) {
              iCurrentNode = kThisLeftNode[iCurrentNode];
            } else {
              iCurrentNode = kThisRightNode[iCurrentNode];
            }
          } else  // categorical
          {
            const Rcpp::IntegerVector kMySplits =
                kSplits[kThisSplitCode[iCurrentNode]];
            if (kMySplits.size() < (int)dX + 1) {
              iCurrentNode = kThisMissingNode[iCurrentNode];
            } else {
              const int iCatSplitIndicator = kMySplits[(int)dX];
              if (iCatSplitIndicator == -1) {
                iCurrentNode = kThisLeftNode[iCurrentNode];
              } else if (iCatSplitIndicator == 1) {
                iCurrentNode = kThisRightNode[iCurrentNode];
              } else  // categorical level not present in training
              {
                iCurrentNode = kThisMissingNode[iCurrentNode];
              }
            }
          }
        }
        predicted_func[kNumCovarRows * prediction_iteration +
                       kNumCovarRows * class_num + obs_num] +=
            kThisSplitCode[iCurrentNode];  // add the prediction
      }                                    // iObs
      tree_num++;
    }  // iTree
  }    // iPredIteration

  return Rcpp::wrap(predicted_func);
  END_RCPP
}

//-----------------------------------
// Function: gbm_plot
//
// Returns: predicted function - a Rcpp::NumericVector
//
// Description: Calculates predictions using only
//				effects of specific variables.
//
// Parameters:
//  covariates - SEXP containing the predictor values - becomes
//				  const Rcpp::NumericMatrix.
//  whichvar- SEXP containing indices specifying which var columns of the
//			  covariates are important - const Rcpp::IntegerVector.
//  num_trees - SEXP containing an int or vector of ints specifying the number
//  of
//				  trees to make predictions on - stored as const
// Rcpp::IntegerVector.
//				  trees to make predictions on - stored as const
// Rcpp::IntegerVector.
//  init_func_est - SEXP specifying initial prediction estimate - double.
//  fitted_trees - SEXP containing lists defining the previously fitted trees -
//					stored as const Rcpp::GenericVector.
//  categorical_splits - SEXP containing  list of the categories of the split
//  variables
//						defining a tree - stored as a
// const
// Rcpp::GenericVector.
//  var_types -  SEXP containing integers specifying whether the variables of
//  interest are
//				continuous/nominal- stored as const
// Rcpp::IntegerVector.
//-----------------------------------

SEXP gbm_plot(
    SEXP covariates,          // vector or matrix of points to make predictions
    SEXP whichvar,            // index of which var cols of X are
    SEXP num_trees,           // number of trees to use
    SEXP init_func_est,       // initial value
    SEXP fitted_trees,        // tree list object
    SEXP categorical_splits,  // categorical split list object
    SEXP var_types            // vector of variable types
    ) {
  BEGIN_RCPP
  int tree_num = 0;
  int obs_num = 0;
  int class_num = 0;
  const Rcpp::NumericMatrix kCovarMat(covariates);
  const int kNumRows = kCovarMat.nrow();
  const int kNumTrees = Rcpp::as<int>(num_trees);
  const Rcpp::IntegerVector kWhichVars(whichvar);
  const Rcpp::GenericVector kFittedTrees(fitted_trees);
  const Rcpp::GenericVector kSplits(categorical_splits);
  const Rcpp::IntegerVector kVarType(var_types);

  Rcpp::NumericVector predicted_func(kNumRows, Rcpp::as<double>(init_func_est));

  if (kCovarMat.ncol() != kWhichVars.size()) {
    throw gbm_exception::InvalidArgument("shape mismatch");
  }

  for (tree_num = 0; tree_num < kNumTrees; tree_num++) {
    const Rcpp::GenericVector kThisTree = kFittedTrees[tree_num];
    const Rcpp::IntegerVector kThisSplitVar = kThisTree[0];
    const Rcpp::NumericVector kThisSplitCode = kThisTree[1];
    const Rcpp::IntegerVector kThisLeftNode = kThisTree[2];
    const Rcpp::IntegerVector kThisRightNode = kThisTree[3];
    const Rcpp::IntegerVector kThisMissingNode = kThisTree[4];
    const Rcpp::NumericVector kThisWeight = kThisTree[6];
    for (obs_num = 0; obs_num < kNumRows; obs_num++) {
      NodeStack stack;
      stack.push(0, 1.0);
      while (!stack.empty()) {
        const std::pair<int, double> top = stack.pop();
        int iCurrentNode = top.first;
        const double kWeight = top.second;

        if (kThisSplitVar[iCurrentNode] == -1)  // terminal node
        {
          predicted_func[class_num * kNumRows + obs_num] +=
              kWeight * kThisSplitCode[iCurrentNode];
        } else  // non-terminal node
        {
          // is this a split variable that interests me?
          const Rcpp::IntegerVector::const_iterator found =
              std::find(kWhichVars.begin(), kWhichVars.end(),
                        kThisSplitVar[iCurrentNode]);

          if (found != kWhichVars.end()) {
            const int kPredVar = found - kWhichVars.begin();
            const double kXValue = kCovarMat(obs_num, kPredVar);
            // missing?
            if (ISNA(kXValue)) {
              stack.push(kThisMissingNode[iCurrentNode], kWeight);
            }
            // continuous?
            else if (kVarType[kThisSplitVar[iCurrentNode]] == 0) {
              if (kXValue < kThisSplitCode[iCurrentNode]) {
                stack.push(kThisLeftNode[iCurrentNode], kWeight);
              } else {
                stack.push(kThisRightNode[iCurrentNode], kWeight);
              }
            } else  // categorical
            {
              const Rcpp::IntegerVector kCatSplits =
                  kSplits[kThisSplitCode[iCurrentNode]];

              const int kCatSplitIndicator = kCatSplits[kXValue];
              if (kCatSplitIndicator == -1) {
                stack.push(kThisLeftNode[iCurrentNode], kWeight);
              } else if (kCatSplitIndicator == 1) {
                stack.push(kThisRightNode[iCurrentNode], kWeight);
              } else  // handle unused level
              {
                stack.push(kThisMissingNode[iCurrentNode], kWeight);
              }
            }
          }     // iPredVar != -1
          else  // not interested in this split, average left and right
          {
            const int kRight = kThisRightNode[iCurrentNode];
            const int kLeft = kThisLeftNode[iCurrentNode];
            const double kRightWeight = kThisWeight[kRight];
            const double kLeftWeight = kThisWeight[kLeft];
            stack.push(kRight,
                       kWeight * kRightWeight / (kRightWeight + kLeftWeight));
            stack.push(kLeft,
                       kWeight * kLeftWeight / (kRightWeight + kLeftWeight));
          }
        }  // non-terminal node
      }    // while(cStackNodes > 0)
    }      // iObs
  }        // iTree

  return Rcpp::wrap(predicted_func);
  END_RCPP
}  // gbm_plot

}  // end extern "C"
