//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//
//  File:       node_search.h
//
//  License:    GNU GPL (version 2 or later)
//
//  Contents:   does the searching for where to split a node
//
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef NODESEARCH_H
#define NODESEARCH_H

//------------------------------
// Includes
//------------------------------
#include "databag.h"
#include "dataset.h"
#include "node.h"
#include "node_parameters.h"
#include "vec_varsplitters.h"
#include "vec_nodeparams.h"
#include "parallel_details.h"

#include <vector>

//------------------------------
// Class Definition
//------------------------------
class CNodeSearch {
 public:
  //----------------------
  // Public Constructors
  //----------------------
  CNodeSearch(unsigned long treedepth, unsigned long minobs,
              const parallel_details& parallel);

  //---------------------
  // Public destructor
  //---------------------
  ~CNodeSearch();

  //---------------------
  // Public Functions
  //---------------------
  void GenerateAllSplits(std::vector<CNode*>& term_nodes_ptrs, 
                         const CDataset& kData,
                         const Bag& kBag, 
                         const std::vector<double>& residuals,
                         std::vector<unsigned long>& data_node_assigns);
  double CalcImprovementAndSplit(std::vector<CNode*>& term_nodes_ptrs,
                                 const CDataset& kData,
                                 std::vector<unsigned long>& data_node_assigns);

  int get_num_threads() const { return parallel_.get_num_threads(); }
  int get_array_chunk_size() const { return parallel_.get_array_chunk_size(); }
  
 private:
  //---------------------
  // Private Functions
  //---------------------
  void ReassignData(unsigned long splittednode_index,
                    std::vector<CNode*>& term_nodes_ptrs, 
                    const CDataset& kData,
                    std::vector<unsigned long>& data_node_assigns);

  //---------------------
  // Private Variables
  //---------------------
  // Best Splits
  VecNodeParams best_splits_;

  // Number of terminal nodes
  unsigned long num_terminal_nodes_;
  unsigned long min_num_node_obs_;

  // parallelization
  parallel_details parallel_;
};

#endif  // NODESEARCH_H
