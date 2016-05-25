//------------------------------------------------------------------------------
//
//  File:       terminalStrategy.h
//
//  Description: Strategies for terminal nodes
//
//	Author: 	James Hickey
//------------------------------------------------------------------------------

#ifndef TERMINALSTRATEGY_H
#define TERMINALSTRATEGY_H

//------------------------------
// Includes
//------------------------------
#include "dataset.h"
#include "generic_node_strategy.h"
#include "node.h"
#include <Rcpp.h>

//------------------------------
// Class Definition
//------------------------------
class TerminalStrategy : public GenericNodeStrategy {
 public:
  //----------------------
  // Public Constructors
  //----------------------
  TerminalStrategy(CNode* node) : nodecontext_(node){};

  //---------------------
  // Public destructor
  //---------------------
  ~TerminalStrategy() { nodecontext_ = NULL; };

  //---------------------
  // Public Functions
  //---------------------
  void Adjust(unsigned long min_num_node_obs) { return; };
  void Predict(const CDataset& kData, unsigned long rownum,
               double& delta_estimate) {
    delta_estimate = nodecontext_->prediction;
  }
  void GetVarRelativeInfluence(double* relative_influence) { return; }
  void PrintSubTree(unsigned long indent) {
    for (unsigned long i = 0; i < indent; i++) Rprintf("  ");
    Rprintf("N=%f, Prediction=%f *\n", nodecontext_->totalweight,
            nodecontext_->prediction);
  }
  signed char WhichNode(const CDataset& kData, unsigned long obs_num) {
    signed char returnvalue = 0;
    double xval = kData.x_value(obs_num, nodecontext_->split_var);
    if (!ISNA(xval)) {
      if (xval < nodecontext_->splitvalue) {
        returnvalue = -1;
      } else {
        returnvalue = 1;
      }
    }
    // if missing value returns 0

    return returnvalue;
  }
  void TransferTreeToRList(int& nodeid, const CDataset& kData, int* splitvar,
                           double* splitpoint, int* leftnodes, int* rightnodes,
                           int* missingnodes, double* error_reduction,
                           double* weights, double* predictions,
                           VecOfVectorCategories& splitcodes_vec,
                           int prev_categorical_splits, double shrinkage) {
    splitvar[nodeid] = -1;
    splitpoint[nodeid] = shrinkage * nodecontext_->prediction;
    leftnodes[nodeid] = -1;
    rightnodes[nodeid] = -1;
    missingnodes[nodeid] = -1;
    error_reduction[nodeid] = 0.0;
    weights[nodeid] = nodecontext_->totalweight;
    predictions[nodeid] = shrinkage * nodecontext_->prediction;

    nodeid++;
  }

 private:
  CNode* nodecontext_;
};
#endif  // TERMINALSTRATEGY_H
