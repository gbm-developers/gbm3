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
#include "node_parameters.h"

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

    inline void reset(){ cTerminalNodes = 1; }
    void set_search_rootnode(CNode& rootNode){ variableSplitters[0].Set(rootNode); }

private:
	//---------------------
	// Private Functions
	//---------------------
    void ReAssignData(unsigned long splittedNodeIndex, vector<CNode*>& vecpTermNodes,
    					const CDataset& data, vector<unsigned long>& aiNodeAssign);

	//---------------------
	// Private Variables
	//---------------------
    // Splitters for variable sets
    std::vector<VarSplitter> variableSplitters;

    // Number of terminal nodes
    unsigned long cTerminalNodes;
    unsigned long minNumObs;
};

#endif // NODESEARCH_H
