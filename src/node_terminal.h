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
    void Adjust(unsigned long cMinObsInNode);

    void PrintSubtree(unsigned long cIndent);
    void TransferTreeToRList(int &iNodeID,
			     const CDataset &data,
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
    
    void ApplyShrinkage(double dLambda);
    void Predict(const CDataset &data,
		 unsigned long i,
		 double &dFadj);
    void Predict(double *adX,
		 unsigned long cRow,
		 unsigned long cCol,
		 unsigned long iRow,
		 double &dFadj);

    void GetVarRelativeInfluence(double *adRelInf);
    void RecycleSelf(CNodeFactory *pNodeFactory);
};

typedef CNodeTerminal *PCNodeTerminal;
typedef vector<PCNodeTerminal> VEC_P_NODETERMINAL;
#endif // NODETERMINAL_H



