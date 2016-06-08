//------------------------------------------------------------------------------
//
//  File:       continuousStrategy.h
//
//  Description: strategies for continuous splits.
//
//	Author: 	James Hickey
//------------------------------------------------------------------------------

#ifndef CONTINUOUSSTRATEGY_H
#define CONTINUOUSSTRATEGY_H

//------------------------------
// Includes
//------------------------------
#include "dataset.h"
#include "node.h"
#include "generic_node_strategy.h"
#include <Rcpp.h>

//------------------------------
// Class Definition
//------------------------------
class ContinuousStrategy : public GenericNodeStrategy {
 public:
  //----------------------
  // Public Constructors
  //----------------------
  ContinuousStrategy(CNode* node) : node_context_(node) { is_split_=true; };

  //---------------------
  // Public destructor
  //---------------------
  ~ContinuousStrategy() { node_context_ = NULL; };

  //---------------------
  // Public Functions
  //---------------------
  void Adjust(unsigned long min_num_node_obs) {
    node_context_->left_child()->Adjust(min_num_node_obs);
    node_context_->right_child()->Adjust(min_num_node_obs);

    if ((node_context_->missing_child()->is_terminal()) &&
        (node_context_->missing_child()->get_numobs() < min_num_node_obs)) {
      node_context_->set_prediction(
          ((node_context_->left_child()->get_totalweight()) *
               (node_context_->left_child()->get_prediction()) +
           (node_context_->right_child()->get_totalweight()) *
               (node_context_->right_child()->get_prediction())) /
          (node_context_->left_child()->get_totalweight() +
           node_context_->right_child()->get_totalweight()));
      node_context_->missing_child()->set_prediction(
          node_context_->get_prediction());
    } else {
      node_context_->missing_child()->Adjust(min_num_node_obs);
      node_context_->set_prediction(
          ((node_context_->left_child()->get_totalweight()) *
               (node_context_->left_child()->get_prediction()) +
           (node_context_->right_child()->get_totalweight()) *
               (node_context_->right_child()->get_prediction()) +
           (node_context_->missing_child()->get_totalweight()) *
               (node_context_->missing_child()->get_prediction())) /
          (node_context_->left_child()->get_totalweight() +
           node_context_->right_child()->get_totalweight() +
           node_context_->missing_child()->get_totalweight()));
    }
  }
  void Predict(const CDataset& kData, unsigned long row_num,
               double& deltafunc_est) {
    signed char whichnode = ContinuousStrategy::WhichNode(kData, row_num);
    if (whichnode == -1) {
      node_context_->left_child()->Predict(kData, row_num, deltafunc_est);
    } else if (whichnode == 1) {
      node_context_->right_child()->Predict(kData, row_num, deltafunc_est);
    } else {
      node_context_->missing_child()->Predict(kData, row_num, deltafunc_est);
    }
  }
  void GetVarRelativeInfluence(double* relative_influence) {
    relative_influence[node_context_->get_split_var()] +=
        node_context_->get_improvement();
    node_context_->left_child()->GetVarRelativeInfluence(relative_influence);
    node_context_->right_child()->GetVarRelativeInfluence(relative_influence);
  }
  void PrintSubTree(unsigned long indent) {
    const std::size_t leftcategory = node_context_->get_leftcategory().size();

    for (unsigned long i = 0; i < indent; i++) Rprintf("  ");
    Rprintf("N=%f, Improvement=%f, Prediction=%f, NA pred=%f\n",
            node_context_->get_totalweight(), node_context_->get_improvement(),
            node_context_->get_prediction(),
            (node_context_->missing_child().get() == 0
                 ? 0.0
                 : node_context_->missing_child()->get_prediction()));

    for (unsigned long i = 0; i < indent; i++) Rprintf("  ");
    Rprintf("V%d in ", node_context_->get_split_var());
    for (unsigned long i = 0; i < leftcategory; i++) {
      Rprintf("%d", node_context_->get_leftcategory()[i]);
      if (i < leftcategory - 1) Rprintf(",");
    }
    Rprintf("\n");
    node_context_->left_child()->PrintSubtree((indent + 1));

    for (unsigned long i = 0; i < indent; i++) Rprintf("  ");
    Rprintf("V%d not in ", node_context_->get_split_var());
    for (unsigned long i = 0; i < leftcategory; i++) {
      Rprintf("%d", node_context_->get_leftcategory()[i]);
      if (i < leftcategory - 1) Rprintf(",");
    }
    Rprintf("\n");
    node_context_->right_child()->PrintSubtree(indent + 1);

    for (unsigned long i = 0; i < indent; i++) Rprintf("  ");
    Rprintf("missing\n");
    node_context_->missing_child()->PrintSubtree(indent + 1);
  }
  signed char WhichNode(const CDataset& kData, unsigned long obs_num) {
    signed char return_value = 0;
    double xval = kData.x_value(obs_num, node_context_->get_split_var());
    if (!ISNA(xval)) {
      if (xval < node_context_->get_splitvalue()) {
        return_value = -1;
      } else {
        return_value = 1;
      }
    }
    // if missing value returns 0

    return return_value;
  }
  void TransferTreeToRList(int& nodeid, const CDataset& kData, int* splitvar,
                           double* splitvalues, int* leftnodes, int* rightnodes,
                           int* missingnodes, double* error_reduction,
                           double* weights, double* predictions,
                           VecOfVectorCategories& splitcodes_vec,
                           int prev_categorical_splits, double shrinkage) {
    int thisnode_id = nodeid;
    splitvar[thisnode_id] = node_context_->get_split_var();
    splitvalues[thisnode_id] = node_context_->get_splitvalue();
    error_reduction[thisnode_id] = node_context_->get_improvement();
    weights[thisnode_id] = node_context_->get_totalweight();
    predictions[thisnode_id] = shrinkage * node_context_->get_prediction();

    nodeid++;
    leftnodes[thisnode_id] = nodeid;
    node_context_->left_child()->TransferTreeToRList(
        nodeid, kData, splitvar, splitvalues, leftnodes, rightnodes,
        missingnodes, error_reduction, weights, predictions, splitcodes_vec,
        prev_categorical_splits, shrinkage);

    rightnodes[thisnode_id] = nodeid;
    node_context_->right_child()->TransferTreeToRList(
        nodeid, kData, splitvar, splitvalues, leftnodes, rightnodes,
        missingnodes, error_reduction, weights, predictions, splitcodes_vec,
        prev_categorical_splits, shrinkage);
    missingnodes[thisnode_id] = nodeid;
    node_context_->missing_child()->TransferTreeToRList(
        nodeid, kData, splitvar, splitvalues, leftnodes, rightnodes,
        missingnodes, error_reduction, weights, predictions, splitcodes_vec,
        prev_categorical_splits, shrinkage);
  }

 private:
  CNode* node_context_;
};
#endif  // CONTINUOUSSTRATEGY_H
