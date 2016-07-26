//  GBM by Greg Ridgeway  Copyright (C) 2003
//-----------------------------------
// Includes
//-----------------------------------
#include <algorithm>
#include "tree.h"

//----------------------------------------
// Function Members - Public
//----------------------------------------
CCARTTree::CCARTTree(const TreeParams& treeconfig)
    : min_num_node_obs_(treeconfig.min_obs_in_node),
      kTreeDepth_(treeconfig.depth),
      kShrinkage_(treeconfig.shrinkage),
      error_(0.0),
      totalnodecount_(1),
      rootnode_(),
      terminalnode_ptrs_(2 * kTreeDepth_ + 1, 0),
      data_node_assignment_(treeconfig.num_trainrows, 0),
      parallel_(treeconfig.parallel) {
  if (kTreeDepth_ < 1) {
    throw gbm_exception::InvalidArgument();
  }
}

//------------------------------------------------------------------------------
// Grows a regression tree
//------------------------------------------------------------------------------
void CCARTTree::Grow(const std::vector<double>& residuals,
		     const CDataset& kData,
                     const Bag& kBag,
                     const std::vector<double>& kDeltaEstimate) {
  if ((residuals.size() < kData.get_trainsize()) ||
      (kDeltaEstimate.size() < kData.get_trainsize())) {
    throw gbm_exception::InvalidArgument();
  }

  double sumz = 0.0;
  double sum_zsquared = 0.0;
  double totalw = 0.0;

  // Move to data -- FOR TIME BEING
  for (unsigned long obs_num = 0; obs_num < kData.get_trainsize(); obs_num++) {
    //if (kBag.get_element(obs_num)) {
      // get the initial sums and sum of squares and total weight
      sumz += kData.weight_ptr()[obs_num] * residuals[obs_num];
      sum_zsquared +=
          kData.weight_ptr()[obs_num] * residuals[obs_num] * residuals[obs_num];
      totalw += kData.weight_ptr()[obs_num];
    //}
  }

  error_ = sum_zsquared - sumz * sumz / totalw;
  rootnode_.reset(new CNode(NodeDef(sumz, totalw, kBag.get_total_in_bag())));
  terminalnode_ptrs_[0] = rootnode_.get();
  CNodeSearch new_node_searcher(kTreeDepth_, min_num_node_obs_, parallel_);

  // build the tree structure
  for (long cDepth = 0; cDepth < kTreeDepth_; cDepth++) {
    // Generate all splits
    new_node_searcher.GenerateAllSplits(terminalnode_ptrs_, kData, kBag,
                                        residuals, data_node_assignment_);
    double bestImprov = new_node_searcher.CalcImprovementAndSplit(
        terminalnode_ptrs_, kData, data_node_assignment_);

    // Make the best split if possible
    if (bestImprov <= 0) {
      break;
    }
    // setup the new nodes and add them to the tree

    totalnodecount_ += 3;

  }  // end tree growing
  // throw gbm_exception::Failure("Here");
  // DEBUG
  // Print();
}

void CCARTTree::PredictValid(const CDataset& kData,
                             unsigned long num_validation_points,
                             std::vector<double>& delta_estimates) {
  unsigned int i = 0;
  for (i = kData.nrow() - num_validation_points; i < kData.nrow(); i++) {
    rootnode_->Predict(kData, i, delta_estimates[i]);
    delta_estimates[i] *= kShrinkage_;
  }
}

void CCARTTree::Adjust(std::vector<double>& delta_estimates) {
  rootnode_->Adjust(min_num_node_obs_);

  // predict for the training observations
  for (unsigned long obs_num = 0; obs_num < data_node_assignment_.size();
       obs_num++) {
    delta_estimates[obs_num] =
        terminalnode_ptrs_[data_node_assignment_[obs_num]]->get_prediction();
  }
}

void CCARTTree::Print() {
  if (rootnode_.get() != 0) {
    rootnode_->PrintSubtree(0);
    Rprintf("shrinkage: %f\n", kShrinkage_);
    Rprintf("initial error: %f\n\n", error_);
  }
}

//-----------------------------------
// Function: TransferTreeToRList
//
// Returns: none
//
// Description:
//
// Parameters:
//
//-----------------------------------
void CCARTTree::TransferTreeToRList(const CDataset& kData, int* splitvar,
                                    double* splitvalues, int* leftnodes,
                                    int* rightnodes, int* missingnodes,
                                    double* error_reduction, double* weights,
                                    double* predictions,
                                    VecOfVectorCategories& splitcodes_vec,
                                    int prev_categorical_splits) {
  int nodeid = 0;
  if (rootnode_.get() != 0) {
    rootnode_->TransferTreeToRList(
        nodeid, kData, splitvar, splitvalues, leftnodes, rightnodes,
        missingnodes, error_reduction, weights, predictions, splitcodes_vec,
        prev_categorical_splits, kShrinkage_);
  } else {
    throw gbm_exception::Failure(
        "Can't transfer to list - RootNode does not exist.");
  }
}
