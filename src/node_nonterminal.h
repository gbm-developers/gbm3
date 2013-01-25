//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//
//  File:       node_nonterminal.h
//
//  License:    GNU GPL (version 2 or later)
//
//  Contents:   a node in the tree
//        	  
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef NODENONTERMINAL_H
#define NODENONTERMINAL_H

#include "node.h"
#include "node_terminal.h"

class CNodeNonterminal : public CNode
{
public:

    CNodeNonterminal();
    virtual ~CNodeNonterminal();
    virtual GBMRESULT Adjust(unsigned long cMinObsInNode);

    virtual signed char WhichNode(CDataset *pData,
                                  unsigned long iObs) = 0;
    virtual signed char WhichNode(double *adX,
                                  unsigned long cRow,
                                  unsigned long cCol,
                                  unsigned long iRow) = 0;
    virtual GBMRESULT TransferTreeToRList(int &iNodeID,
                                        CDataset *pData,
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
                                        double dShrinkage) = 0;

    GBMRESULT Predict(CDataset *pData, 
                    unsigned long iRow, 
                    double &dFadj);
    GBMRESULT Predict(double *adX,
                    unsigned long cRow,
                    unsigned long cCol,
                    unsigned long iRow,
                    double &dFadj);

    GBMRESULT GetVarRelativeInfluence(double *adRelInf);
    virtual GBMRESULT RecycleSelf(CNodeFactory *pNodeFactory) = 0;

    CNode *pLeftNode;
    CNode *pRightNode;
    CNode *pMissingNode;
    unsigned long iSplitVar;
    double dImprovement;
};

typedef CNodeNonterminal *PCNodeNonterminal;

#endif // NODENONTERMINAL_H



