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
#include "node_continuous.h"
#include "node_categorical.h"

using namespace std;

class CNodeSearch
{
public:

    CNodeSearch();
    ~CNodeSearch();
    void Initialize(unsigned long cMinObsInNode);

    void IncorporateObs(double dX,
			double dZ,
			double dW,
			long lMonotone);

    void Set(double dSumZ,
	     double dTotalW,
	     unsigned long cTotalN,
	     CNode* pThisNode,
	     CNode **ppParentPointerToThisNode);
    void ResetForNewVar(unsigned long iWhichVar,
			long cVarClasses);
    
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

    double BestImprovement() { return dBestImprovement; }
    void SetToSplit()
    {
        fIsSplit = true;
    };
    void SetupNewNodes(CNode* &pNewSplitNode,
		       CNode* &pNewLeftNode,
		       CNode* &pNewRightNode,
		       CNode* &pNewMissingNode);

    void EvaluateCategoricalSplit();
    void WrapUpCurrentVariable();
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

    std::vector<double> adGroupSumZ;
    std::vector<double> adGroupW;
    std::vector<unsigned long> acGroupN;
    std::vector<double> adGroupMean;
    // this is an int to fit in with R API
    // it's probably best not to ask.
    std::vector<int> aiCurrentCategory;
    std::vector<unsigned long> aiBestCategory;

    CNode* pThisNode;
    CNode **ppParentPointerToThisNode;
};

#endif // NODESEARCH_H
