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
    CCARTTree(double shrinkage=1.0, long depth = 1);

	//---------------------
	// Public destructor
	//---------------------
    ~CCARTTree();

	//---------------------
	// Public Functions
	//---------------------
    void grow(double *adZ,
	      const CDataset& data,
	      const double *adF,
	      unsigned long cMinObsInNode,
	      std::vector<unsigned long>& aiNodeAssign,
	      CNodeSearch& aNodeSearch);
    void Reset();

    CNode* GetRootNode();
    const CNode* GetRootNode() const;

    void PredictValid(const CDataset &pData,
		      unsigned long nValid,
		      double *adFadj);
    void Adjust(const std::vector<unsigned long>& aiNodeAssign,
		double *adFadj,
		unsigned long cMinObsInNode);

    const unsigned long& GetNodeCount() const { return totalnodecount_; }
    vector<CNode*>& GetTermNodes() { return terminalnode_ptrs_; }
    const double& GetShrinkageConst() const { return kShrinkage_; }
    void Print();

private:
	//---------------------
	// Private Variables
	//---------------------
    CNode* rootnode_;
    vector<CNode*> terminalnode_ptrs_;

    const long kTreeDepth_;
    const double kShrinkage_;
    double error_; // total squared error before carrying out the splits
    unsigned long totalnodecount_;

};

#endif // TREE_H
