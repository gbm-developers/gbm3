//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//
//  File:       node_terminal.h
//
//  License:    GNU GPL (version 2 or later)
//
//  Contents:   terminal node class
//        	  
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef NODETERMINAL_H
#define NODETERMINAL_H

#include <vector>
#include "dataset.h"
#include "node.h"

using namespace std;

class CNodeTerminal : public CNode
{
public:

    CNodeTerminal();
    ~CNodeTerminal();
    GBMRESULT Adjust(unsigned long cMinObsInNode);

    GBMRESULT PrintSubtree(unsigned long cIndent);
    GBMRESULT TransferTreeToRList(int &iNodeID,
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
                                double dShrinkage);

    GBMRESULT ApplyShrinkage(double dLambda);
    GBMRESULT Predict(CDataset *pData, 
                    unsigned long i, 
                    double &dFadj);
    GBMRESULT Predict(double *adX,
                    unsigned long cRow,
                    unsigned long cCol,
                    unsigned long iRow,
                    double &dFadj);

    GBMRESULT GetVarRelativeInfluence(double *adRelInf);
    GBMRESULT RecycleSelf(CNodeFactory *pNodeFactory);
};

typedef CNodeTerminal *PCNodeTerminal;
typedef vector<PCNodeTerminal> VEC_P_NODETERMINAL;
#endif // NODETERMINAL_H



