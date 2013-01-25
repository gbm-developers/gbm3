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

#ifndef TREGBM_H
#define TREGBM_H

#include <cstdio>
#include <cfloat>
#include <algorithm>
#include "dataset.h"
#include "node_factory.h"
#include "node_search.h"


class CCARTTree
{
public:

    CCARTTree();
    ~CCARTTree();

    GBMRESULT Initialize(CNodeFactory *pNodeFactory);
    GBMRESULT grow(double *adZ, 
                 CDataset *pData, 
                 double *adAlgW,
                 double *adF,
                 unsigned long nTrain,
                 unsigned long nBagged,
                 double dLambda,
                 unsigned long cMaxDepth,
                 unsigned long cMinObsInNode,
                 bool *afInBag,
                 unsigned long *aiNodeAssign,
                 CNodeSearch *aNodeSearch,
                 VEC_P_NODETERMINAL &vecpTermNodes);
    GBMRESULT Reset();

    GBMRESULT TransferTreeToRList(CDataset *pData,
                                int *aiSplitVar,
                                double *adSplitPoint,
                                int *aiLeftNode,
                                int *aiRightNode,
                                int *aiMissingNode,
                                double *adErrorReduction,
                                double *adWeight,
                                double *adPred,
                                VEC_VEC_CATEGORIES &vecSplitCodes,
                                int cCatSplitsOld,
                                double dShrinkage);

    GBMRESULT PredictValid(CDataset *pData, 
                         unsigned long nValid, 
                         double *adFadj);

    GBMRESULT Predict(double *adX,
                    unsigned long cRow, 
                    unsigned long cCol, 
                    unsigned long iRow, 
                    double &dFadj);
    GBMRESULT Adjust(unsigned long *aiNodeAssign,
                   double *adFadj,
                   unsigned long cTrain,
                   VEC_P_NODETERMINAL &vecpTermNodes,
                   unsigned long cMinObsInNode);

    GBMRESULT GetNodeCount(int &cNodes);
    GBMRESULT SetShrinkage(double dShrink)
    {
        this->dShrink = dShrink;
        return GBM_OK;
    }
    double GetShrinkage() {return dShrink;}

    GBMRESULT Print();
    GBMRESULT GetVarRelativeInfluence(double *adRelInf);


    double dError; // total squared error before carrying out the splits
private:
    GBMRESULT GetBestSplit(CDataset *pData,
                         unsigned long nTrain,
                         CNodeSearch *aNodeSearch,
                         unsigned long cTerminalNodes,
                         unsigned long *aiNodeAssign,
                         bool *afInBag,
                         double *adZ,
                         double *adW,
                         unsigned long &iBestNode,
                         double &dBestNodeImprovement);

    CNode *pRootNode;
    double dShrink;

    // objects used repeatedly
    unsigned long cDepth;
    unsigned long cTerminalNodes;
    unsigned long cTotalNodeCount;
    unsigned long iObs;
    unsigned long iWhichNode;

    unsigned long iBestNode;
    double dBestNodeImprovement;

    double dSumZ;
    double dSumZ2;
    double dTotalW;
    signed char schWhichNode;

    CNodeFactory *pNodeFactory;
    CNodeNonterminal *pNewSplitNode;
    CNodeTerminal *pNewLeftNode;
    CNodeTerminal *pNewRightNode;
    CNodeTerminal *pNewMissingNode;
    CNodeTerminal *pInitialRootNode;
};

typedef CCARTTree *PCCARTTree;


#endif // TREGBM_H



