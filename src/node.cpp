//  GBM by Greg Ridgeway  Copyright (C) 2003

//-----------------------------------
// Includes
//-----------------------------------
#include "node.h"
#include "terminal_strategy.h"
#include "continuous_strategy.h"
#include "categorical_strategy.h"

//----------------------------------------
// Function Members - Public
//----------------------------------------
CNode::CNode(const NodeDef& kDefn)
    : node_strategy_(new TerminalStrategy(this)),
      left_node_ptr_(NULL),
      right_node_ptr_(NULL),
      missing_node_ptr_(NULL),
      split_var_(0.0), improvement_(0.0),
      prediction_(kDefn.prediction()),
      totalweight_(kDefn.get_totalweight()),
      numobs_(kDefn.get_num_obs()),
      leftcategory_(), splitvalue_(0.0), splitdetermined_(false) {}

void CNode::SetStrategy(bool is_continuous_split) {

  if(is_continuous_split) {
	  node_strategy_.reset(new ContinuousStrategy(this));
  } else {
	  node_strategy_.reset(new CategoricalStrategy(this));
  }

}

void CNode::Adjust(unsigned long min_num_node_obs) {
  node_strategy_->Adjust(min_num_node_obs);
}

void CNode::Predict(const CDataset& kData, unsigned long rownum,
                    double& delta_estimate) {
  node_strategy_->Predict(kData, rownum, delta_estimate);
}

void CNode::GetVarRelativeInfluence(double* relative_influence) {
  node_strategy_->GetVarRelativeInfluence(relative_influence);
}

void CNode::PrintSubtree(unsigned long indent) {
  node_strategy_->PrintSubTree(indent);
}

void CNode::SplitNode(const NodeParams& childrenparams) {
  // set up a continuous split
  if (childrenparams.split_class() == 0) {
    SetStrategy(true);
  } else {
    SetStrategy(false);
    // the types are confused here
    leftcategory_.resize(1 + (unsigned long)childrenparams.split_value());
    std::copy(childrenparams.get_ordering().begin(),
              childrenparams.get_ordering().begin() + leftcategory_.size(),
              leftcategory_.begin());
  }

  split_var_ = childrenparams.split_variable();
  splitvalue_ = childrenparams.split_value();
  improvement_ = childrenparams.get_improvement();

  // Check that our nodes are defined
  if(!childrenparams.nodes_have_obs()) {
	  throw gbm_exception::Failure("Best split has no observations!");
  }

  left_node_ptr_.reset(new CNode(childrenparams.get_left_def()));
  right_node_ptr_.reset(new CNode(childrenparams.get_right_def()));
  missing_node_ptr_.reset(new CNode(childrenparams.get_missing_def()));
}

signed char CNode::WhichNode(const CDataset& kData, unsigned long obs_num) {
  return node_strategy_->WhichNode(kData, obs_num);
}

bool CNode::is_terminal() const {
	  return node_strategy_->is_split();
}

void CNode::TransferTreeToRList(int& node_id, const CDataset& kData,
                                int* splivar, double* splitvalues,
                                int* leftnodes, int* rightnodes,
                                int* missingnodes, double* error_reduction,
                                double* weights, double* predictions,
                                VecOfVectorCategories& splitcodes_vec,
                                int prev_categorical_split, double shrinkage) {
  node_strategy_->TransferTreeToRList(
      node_id, kData, splivar, splitvalues, leftnodes, rightnodes, missingnodes,
      error_reduction, weights, predictions, splitcodes_vec,
      prev_categorical_split, shrinkage);
}
