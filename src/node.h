//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//
//  File:       node.h
//
//  License:    GNU GPL (version 2 or later)
//
//  Contents:   a node in the tree
//
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//			   16/03/2016   James Hickey: updated to remove terminal and
//non-terminal nodes
//			   30/03/2016   James Hickey: state pattern to deal with
//continuous and categorical splits.
//------------------------------------------------------------------------------

#ifndef NODE_H
#define NODE_H
//------------------------------
// Includes
//------------------------------
#include "dataset.h"
#include "node_parameters.h"
#include <vector>

//------------------------------
// Class Forwards and Type defs
//------------------------------
class GenericNodeStrategy;

using namespace std;
typedef vector<int> VectorCategories;
typedef vector<VectorCategories> VecOfVectorCategories;

//------------------------------
// Class definition
//------------------------------
class CNode {
 public:
  //----------------------
  // Public Constructors
  //----------------------
  CNode(const NodeDef& kDefn);

  //---------------------
  // Public destructor
  //---------------------
  ~CNode() {};

  //---------------------
  // Public Functions
  //---------------------
  void Adjust(unsigned long min_num_node_obs);
  void Predict(const CDataset& kData, unsigned long rownum,
               double& delta_estimate);

  void GetVarRelativeInfluence(double* relative_influence);
  void SplitNode(NodeParams& childrenparams);
  void PrintSubtree(unsigned long indent);
  void TransferTreeToRList(int& node_iD, const CDataset& kData, int* splitvar,
                           double* splitvalues, int* leftnodes, int* rightnodes,
                           int* missingnodes, double* error_reduction,
                           double* weights, double* predictions,
                           VecOfVectorCategories& splitcodes_vec,
                           int prev_categorical_splits, double shrinkage);
  signed char WhichNode(const CDataset& kData, unsigned long obs_num);
  std::auto_ptr<CNode>& left_child() {
	  return left_node_ptr_;
  }
  const std::auto_ptr<CNode>& left_child() const {
	  return left_node_ptr_;
  }
  std::auto_ptr<CNode>& right_child() {
	  return right_node_ptr_;
  }
  const std::auto_ptr<CNode>& right_child() const {
	  return right_node_ptr_;
  }
  std::auto_ptr<CNode>& missing_child() {
	  return missing_node_ptr_;
  }
  const std::auto_ptr<CNode>& missing_child() const {
	  return missing_node_ptr_;
  }
  unsigned long get_split_var() const {
	  return split_var_;
  }
  double get_improvement() const {
	  return improvement_;
  }
  double get_splitvalue() const {
	  return splitvalue_;
  }
  double get_prediction() const {
	  return prediction_;
  }
  void set_prediction(double pred_val) {
	  prediction_ = pred_val;
  }
  double get_totalweight() const {
	  return totalweight_;
  }
  unsigned long get_numobs() const {
	  return numobs_;
  }
  std::vector<unsigned long>& get_leftcategory() {
	  return leftcategory_;
  }
  bool is_terminal() const;

 private:
  //---------------------
  // Private Functions
  //---------------------
  void SetStrategy(bool is_continuous_split);

  //---------------------
  // Private Variables
  //---------------------
  std::auto_ptr<GenericNodeStrategy> node_strategy_;

  // Pointers to the Node's children
  std::auto_ptr<CNode> left_node_ptr_;
  std::auto_ptr<CNode> right_node_ptr_;
  std::auto_ptr<CNode> missing_node_ptr_;

  // TODO: Currently most useful in printing out tree
  // This nodes parameters
  unsigned long split_var_;
  double improvement_;

  // Properties defining the node
  double prediction_;
  double totalweight_;    // total training weight in node
  unsigned long numobs_;  // number of training observations in node

  // VARIABLES USED IN NODE SPLITTING
  std::vector<unsigned long> leftcategory_;
  double splitvalue_;
};

#endif  // NODE_H
