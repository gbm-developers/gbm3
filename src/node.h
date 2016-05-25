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
// Class Forwards and Enums
//------------------------------
class GenericNodeStrategy;
enum SplitType { kCategorical, kContinuous, kNone };

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
  virtual ~CNode();

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

  //---------------------
  // Public Variables
  //---------------------
  // Pointers to the Node's children
  CNode* left_node_ptr;
  CNode* right_node_ptr;
  CNode* missing_node_ptr;

  // TODO: Currently most useful in printing out tree
  // This nodes parameters
  unsigned long split_var;
  double improvement;

  // Properties defining the node
  double prediction;
  double totalweight;    // total training weight in node
  unsigned long numobs;  // number of training observations in node

  // ENUM FOR strategy
  SplitType splittype;

  // VARIABLES USED IN NODE SPLITTING
  std::vector<unsigned long> leftcategory;
  double splitvalue;

 private:
  //---------------------
  // Private Functions
  //---------------------
  void SetStrategy();

  //---------------------
  // Private Variables
  //---------------------
  GenericNodeStrategy* node_strategy_;
};

#endif  // NODE_H
