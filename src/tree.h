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

#ifndef __tree_h__
#define __tree_h__

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

    long GetNodeCount();
    const long GetNodeCount() const;

    vector<CNode*> GetTermNodes(){return vecpTermNodes;}
    const double GetShrinkageConst() const { return shrinkageConst;}
    void Print();

private:
	//---------------------
	// Private Variables
	//---------------------
    CNode* pRootNode;
    vector<CNode*> vecpTermNodes;

    const long depthOfTree;
    const double shrinkageConst;
    double dError; // total squared error before carrying out the splits
    long cTotalNodeCount;

};

#endif // __tree_h__



