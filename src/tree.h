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
#include <cstdio>
#include <cfloat>
#include <algorithm>
#include <vector>
#include "dataset.h"
#include "node_search.h"
#include <ctime>

//------------------------------
// Class definition
//------------------------------
class CCARTTree
{
public:
	//----------------------
	// Public Constructors
	//----------------------
    CCARTTree(TreeParams treeconfig);

	//---------------------
	// Public destructor
	//---------------------
    ~CCARTTree();

	//---------------------
	// Public Functions
	//---------------------
    void Grow(double* residuals,
	      const CDataset& kData,
	      const double* kFuncEstimate);
    void Reset();

    void PredictValid(const CDataset &kData,
		      unsigned long num_validation_points,
		      double* delta_estimates);
    void Adjust(double* delta_estimates);

    void TransferTreeToRList(const CDataset &kData,
		     int* splitvar,
		     double* splitvalues,
		     int* leftnodes,
		     int* rightnodes,
		     int* missingnodes,
		     double* error_reduction,
		     double* weights,
		     double* predictions,
		     VEC_VEC_CATEGORIES &splitcodes_vec,
		     int prev_categorical_splits);
    void Print();

    std::vector<unsigned long>& get_node_assignments()
	{
		return data_node_assignment_;
	}
	vector<CNode*>& get_terminal_nodes()
	{
		return terminalnode_ptrs_;
	}
	const double& get_shrinkage_factor() const
	{
		return kShrinkage_;
	}
	const unsigned long& min_num_obs_required() const
	{
		return min_num_node_obs_;
	}
	const unsigned long& size_of_tree() const
	{
		return totalnodecount_;
	}
	 CNode* get_rootnode()
	{
		return rootnode_;
	}
	const CNode* get_rootnode() const
	{
		return rootnode_;
	}

private:
	//---------------------
	// Private Variables
	//---------------------
    CNode* rootnode_;
    vector<CNode*> terminalnode_ptrs_;
    vector<unsigned long> data_node_assignment_;
    CNodeSearch new_node_searcher_;

    unsigned long min_num_node_obs_;
    const long kTreeDepth_;
    const double kShrinkage_;
    double error_; // total squared error before carrying out the splits
    unsigned long totalnodecount_;

};

#endif // TREE_H
