//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//
//  File:       node_terminal.cpp
//
//------------------------------------------------------------------------------

#include "node_terminal.h"
#include "node_factory.h"

CNodeTerminal::CNodeTerminal()
{
    isTerminal = true;
}


CNodeTerminal::~CNodeTerminal()
{
    #ifdef NOISY_DEBUG
    Rprintf("terminal destructor\n");
    #endif
}

void CNodeTerminal::Adjust
(
    unsigned long cMinObsInNode
)
{
}

void CNodeTerminal::ApplyShrinkage
(
    double dLambda
)
{
    dPrediction *= dLambda;
}



void CNodeTerminal::Predict
(
    const CDataset &data,
    unsigned long iRow,
    double &dFadj
)
{
    dFadj = dPrediction;
}




void CNodeTerminal::Predict
(
 double *adX,
 unsigned long cRow,
 unsigned long cCol,
 unsigned long iRow,
 double &dFadj
 )
{
  dFadj = dPrediction;
}



void CNodeTerminal::PrintSubtree
(
 unsigned long cIndent
)
{
  unsigned long i = 0;
  
  for(i=0; i< cIndent; i++) Rprintf("  ");
  Rprintf("N=%f, Prediction=%f *\n",
	  dTrainW,
	  dPrediction);
}


void CNodeTerminal::GetVarRelativeInfluence
(
    double *adRelInf
)
{
}


void CNodeTerminal::RecycleSelf(CNodeFactory *pNodeFactory) {
  pNodeFactory->RecycleNode(this);
};


void CNodeTerminal::TransferTreeToRList
(
    int &iNodeID,
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
    double dShrinkage
)
{
  aiSplitVar[iNodeID] = -1;
  adSplitPoint[iNodeID] = dShrinkage*dPrediction;
  aiLeftNode[iNodeID] = -1;
  aiRightNode[iNodeID] = -1;
  aiMissingNode[iNodeID] = -1;
  adErrorReduction[iNodeID] = 0.0;
  adWeight[iNodeID] = dTrainW;
  adPred[iNodeID] = dShrinkage*dPrediction;
  
  iNodeID++;
}


