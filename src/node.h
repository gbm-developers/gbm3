//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//
//  File:       node.h
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

#ifndef NODGBM_H
#define NODGBM_H

#include <vector>
#include "dataset.h"
#include "buildinfo.h"


class CNodeFactory;

using namespace std;

typedef vector<int> VEC_CATEGORIES;
typedef vector<VEC_CATEGORIES> VEC_VEC_CATEGORIES;


class CNode
{
public:

    CNode();
    virtual ~CNode();
    virtual void Adjust(unsigned long cMinObsInNode) = 0;
    virtual void Predict(CDataset &data,
			 unsigned long iRow,
			 double &dFadj) = 0;
    virtual void Predict(double *adX,
			 unsigned long cRow,
			 unsigned long cCol,
			 unsigned long iRow,
			 double &dFadj) = 0;
    static double Improvement
    (
        double dLeftW,
        double dRightW,
        double dMissingW,
        double dLeftSum,
        double dRightSum,
        double dMissingSum
    )
    {
        double dTemp = 0.0;
        double dResult = 0.0;

        if(dMissingW == 0.0)
        {
            dTemp = dLeftSum/dLeftW - dRightSum/dRightW;
            dResult = dLeftW*dRightW*dTemp*dTemp/(dLeftW+dRightW);
        }
        else
        {
            dTemp = dLeftSum/dLeftW - dRightSum/dRightW;
            dResult += dLeftW*dRightW*dTemp*dTemp;
            dTemp = dLeftSum/dLeftW - dMissingSum/dMissingW;
            dResult += dLeftW*dMissingW*dTemp*dTemp;
            dTemp = dRightSum/dRightW - dMissingSum/dMissingW;
            dResult += dRightW*dMissingW*dTemp*dTemp;
            dResult /= (dLeftW + dRightW + dMissingW);
        }

        return dResult;
    }


    virtual void PrintSubtree(unsigned long cIndent) = 0;
    virtual void TransferTreeToRList(int &iNodeID,
				     CDataset &data,
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

    virtual void GetVarRelativeInfluence(double *adRelInf) = 0;
    virtual void RecycleSelf(CNodeFactory *pNodeFactory) = 0;
    
    virtual void reset() { dPrediction = 0; }
    double dPrediction;
    double dTrainW;   // total training weight in node
    unsigned long cN; // number of training observations in node
    bool isTerminal;
};

typedef CNode *PCNode;

#endif // NODGBM_H



