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
#include <vector>
#include "dataset.h"
#include "node.h"
#include "varsplitter.h"
#include "nodeParameters.h"

using namespace std;

//------------------------------
// Class Definition
//------------------------------
class CNodeSearch
{
public:
	//----------------------
	// Public Constructors
	//----------------------
    CNodeSearch(int treeDepth, int numColData, unsigned long minObs);

	//---------------------
	// Public destructor
	//---------------------
    ~CNodeSearch();

	//---------------------
	// Public Functions
	//---------------------
    void GenerateAllSplits(vector<CNode*>& vecpTermNodes, const CDataset& data,
    						double* residuals, vector<unsigned long>& aiNodeAssign);
    double CalcImprovementAndSplit(vector<CNode*>& vecpTermNodes, const CDataset& data,
    		vector<unsigned long>& aiNodeAssign);

    inline void Reset(){ cTerminalNodes = 1; }
    void SetRootNode(CNode& rootNode){ variableSplitters[0].Set(rootNode); }

private:
	//---------------------
	// Private Functions
	//---------------------
    void ReAssignData(long splittedNodeIndex, vector<CNode*>& vecpTermNodes,
    					const CDataset& data, vector<unsigned long>& aiNodeAssign);

	//---------------------
	// Private Variables
	//---------------------
    // Splitters for variable sets
    std::vector<VarSplitter> variableSplitters;

    // Number of terminal nodes
    long cTerminalNodes;
    unsigned long minNumObs;
    long totalCache;
};

#endif // NODESEARCH_H

