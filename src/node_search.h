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

#include <vector>
#include "dataset.h"
#include "node.h"
#include "varsplitter.h"
#include "nodeParameters.h"


using namespace std;

class CNodeSearch
{
public:

    CNodeSearch(int numColData, unsigned long minObs);
    ~CNodeSearch();
    void GenerateAllSplits(vector<CNode*>& vecpTermNodes, const CDataset& data,
    						double* residuals, vector<unsigned long>& aiNodeAssign);
    double SplitAndCalcImprovement(vector<CNode*>& vecpTermNodes,
    					const CDataset& data,
    					vector<unsigned long>& aiNodeAssign);
    void Reset();


private:
    //Private methods
    void ReAssignData(long splittedNodeIndex, vector<CNode*>& vecpTermNodes,
    					const CDataset& data, vector<unsigned long>& aiNodeAssign);
    void AssignToNode(CNode& terminalNode);

    // Split Parameters -
    std::vector<VarSplitter> variableSplitters;

    // Number of terminal nodes
    long cTerminalNodes;
};

#endif // NODESEARCH_H
