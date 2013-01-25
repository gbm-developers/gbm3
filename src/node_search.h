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

#include "node_factory.h"
#include "dataset.h"

using namespace std;

class CNodeSearch
{
public:

    CNodeSearch();
    ~CNodeSearch();
    GBMRESULT Initialize(unsigned long cMinObsInNode);

    GBMRESULT IncorporateObs(double dX,
                             double dZ,
                             double dW,
                             long lMonotone);

    GBMRESULT Set(double dSumZ,
                double dTotalW,
                unsigned long cTotalN,
                CNodeTerminal *pThisNode,
                CNode **ppParentPointerToThisNode,
                CNodeFactory *pNodeFactory);
    GBMRESULT ResetForNewVar(unsigned long iWhichVar,
                           long cVarClasses);

    double BestImprovement() { return dBestImprovement; }
    GBMRESULT SetToSplit() 
    {    
        fIsSplit = true;
        return GBM_OK;
    };
    GBMRESULT SetupNewNodes(PCNodeNonterminal &pNewSplitNode,
                          PCNodeTerminal &pNewLeftNode,
                          PCNodeTerminal &pNewRightNode,
                          PCNodeTerminal &pNewMissingNode);

    GBMRESULT EvaluateCategoricalSplit();
    GBMRESULT WrapUpCurrentVariable();
    double ThisNodePrediction() {return pThisNode->dPrediction;}
    bool operator<(const CNodeSearch &ns) {return dBestImprovement<ns.dBestImprovement;}

    unsigned long iBestSplitVar;
    double dBestSplitValue;

    double dBestLeftSumZ;
    double dBestLeftTotalW;
    unsigned long cBestLeftN;

    double dBestRightSumZ;
    double dBestRightTotalW;
    unsigned long cBestRightN;

    double dBestMissingSumZ;
    double dBestMissingTotalW;
    unsigned long cBestMissingN;

    double dCurrentMissingSumZ;
    double dCurrentMissingTotalW;
    unsigned long cCurrentMissingN;

    long cCurrentVarClasses;

    unsigned long iRank;
    double dInitTotalW;
    double dInitSumZ;
    unsigned long cInitN;
    double dBestImprovement;

private:
    bool fIsSplit;

    unsigned long cMinObsInNode;

    long cBestVarClasses;

    double dCurrentLeftSumZ;
    double dCurrentLeftTotalW;
    unsigned long cCurrentLeftN;
    double dCurrentRightSumZ;
    double dCurrentRightTotalW;
    unsigned long cCurrentRightN;
    double dCurrentImprovement;
    unsigned long iCurrentSplitVar;
    double dCurrentSplitValue;

    double dLastXValue;

    double *adGroupSumZ;
    double *adGroupW;
    unsigned long *acGroupN;
    double *adGroupMean;
    int *aiCurrentCategory;
    unsigned long *aiBestCategory;
    const unsigned long k_cMaxClasses;

    CNodeTerminal *pThisNode;
    CNode **ppParentPointerToThisNode;
    CNodeFactory *pNodeFactory;
};

typedef CNodeSearch *PCNodeSearch;

#endif // NODESEARCH_H
