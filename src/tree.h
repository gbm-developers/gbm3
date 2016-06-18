//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//
//  File:       tree.h
//
//  License:    GNU GPL (version 2 or later)
//
//  Contents:   regression tree
//
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef TREE_H
#define TREE_H

//------------------------------
// Includes
//------------------------------
#include "databag.h"
#include "dataset.h"
#include "node_search.h"
#include "parallel_details.h"
#include "treeparams.h"
#include <algorithm>
#include <ctime>
#include <cstdio>
#include <cfloat>
#include <vector>

//------------------------------
// Class definition
//------------------------------
class CCARTTree {
 public:
  //----------------------
  // Public Constructors
  //----------------------
  CCARTTree(const TreeParams& treeconfig);

  //---------------------
  // Public destructor
  //---------------------
  ~CCARTTree(){};

  //---------------------
  // Public Functions
  //---------------------
  void Grow(std::vector<double>& residuals, const CDataset& kData,
            const Bag& kBag, const std::vector<double>& kDeltaEstimate);

  void PredictValid(const CDataset& kData, unsigned long num_validation_points,
                    std::vector<double>& delta_estimates);
  void Adjust(std::vector<double>& delta_estimates);

  void TransferTreeToRList(const CDataset& kData, int* splitvar,
                           double* splitvalues, int* leftnodes, int* rightnodes,
                           int* missingnodes, double* error_reduction,
                           double* weights, double* predictions,
                           VecOfVectorCategories& splitcodes_vec,
                           int prev_categorical_splits);
  void Print();

  std::vector<unsigned long>& get_node_assignments() {
    return data_node_assignment_;
  }
  vector<CNode*>& get_terminal_nodes() { return terminalnode_ptrs_; }
  const double& get_shrinkage_factor() const { return kShrinkage_; }
  const unsigned long& min_num_obs_required() const {
    return min_num_node_obs_;
  }
  const unsigned long& size_of_tree() const { return totalnodecount_; }

  int get_num_threads() const { return parallel_.get_num_threads(); }

 private:
  //---------------------
  // Private Variables
  //---------------------
  unsigned long min_num_node_obs_;
  const long kTreeDepth_;
  const double kShrinkage_;
  double error_;  // total squared error before carrying out the splits
  unsigned long totalnodecount_;

  auto_ptr<CNode> rootnode_;
  vector<CNode*> terminalnode_ptrs_;
  vector<unsigned long> data_node_assignment_;

  parallel_details parallel_;
};

#endif  // TREE_H
