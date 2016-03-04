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
#include <vector>
#include "dataset.h"
#include "node_factory.h"
#include "node_search.h"
#include <ctime>


class CCARTTree
{
public:

    CCARTTree();
    ~CCARTTree();

    void Initialize(CNodeFactory *pNodeFactory);
    void grow(double *adZ,
	      const CDataset &pData,
	      const double *adAlgW,
	      const double *adF,
	      unsigned long nTrain,
	      unsigned long nFeatures,
	      unsigned long nBagged,
	      double dLambda,
	      unsigned long cMaxDepth,
	      unsigned long cMinObsInNode,
	      const bag& afInBag,
	      std::vector<unsigned long>& aiNodeAssign,
	      CNodeSearch *aNodeSearch,
	      VEC_P_NODETERMINAL& vecpTermNodes);
    void Reset();

    void TransferTreeToRList(const CDataset &pData,
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

    void PredictValid(const CDataset &pData,
		      unsigned long nValid,
		      double *adFadj);
    
    void Predict(double *adX,
		 unsigned long cRow,
		 unsigned long cCol,
		 unsigned long iRow,
		 double &dFadj);
    void Adjust(const std::vector<unsigned long>& aiNodeAssign,
		double *adFadj,
		unsigned long cTrain,
		const VEC_P_NODETERMINAL &vecpTermNodes,
		unsigned long cMinObsInNode);
    
    void GetNodeCount(int &cNodes);
    void SetShrinkage(double dShrink)
    {
        this->dShrink = dShrink;
    }
    double GetShrinkage() {return dShrink;}

    void Print();
    void GetVarRelativeInfluence(double *adRelInf);


    double dError; // total squared error before carrying out the splits
private:
    void GetBestSplit(const CDataset &pData,
		      unsigned long nTrain,
		      unsigned long nFeatures,
		      CNodeSearch *aNodeSearch,
		      unsigned long cTerminalNodes,
		      std::vector<unsigned long>& aiNodeAssign,
		      const bag& afInBag,
		      double *adZ,
		      const double *adW,
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



