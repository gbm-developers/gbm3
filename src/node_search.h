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
#include "vec_varsplitters.h"
#include "node_parameters.h"
#include <vector>

using namespace std;

//------------------------------
// Class Definition
//------------------------------
class CNodeSearch {
 public:
  //----------------------
  // Public Constructors
  //----------------------
  CNodeSearch(unsigned long treedepth, unsigned long minobs);

  //---------------------
  // Public destructor
  //---------------------
  ~CNodeSearch();

  //---------------------
  // Public Functions
  //---------------------
  void GenerateAllSplits(vector<CNode* >& term_nodes_ptrs, const CDataset& kData,
		  	  	  	  	 const Bag& kBag,
                         double* residuals,
                         vector<unsigned long>& data_node_assigns);
  double CalcImprovementAndSplit(vector<CNode*>& term_nodes_ptrs,
                                 const CDataset& kData,
                                 vector<unsigned long>& data_node_assigns);

 private:
  //---------------------
  // Private Functions
  //---------------------
  void ReassignData(unsigned long splittednode_index,
                    vector<CNode*>& term_nodes_ptrs, const CDataset& kData,
                    vector<unsigned long>& data_node_assigns);

  //---------------------
  // Private Variables
  //---------------------
  // Splitters for variable sets
   VecVarSplitters best_splits_;

  // Number of terminal nodes
  unsigned long num_terminal_nodes_;
  unsigned long min_num_node_obs_;
};

#endif  // NODESEARCH_H
