//  GBM by Greg Ridgeway  Copyright (C) 2003

#include "node_continuous.h"
#include "node_factory.h"

CNodeContinuous::CNodeContinuous()
{
    dSplitValue = 0.0;
}


CNodeContinuous::~CNodeContinuous()
{
    #ifdef NOISY_DEBUG
    Rprintf("continuous destructor\n");
    #endif
}



void CNodeContinuous::PrintSubtree
(
 unsigned long cIndent
)
{
  unsigned long i = 0;
  
  for(i=0; i< cIndent; i++) Rprintf("  ");
  Rprintf("N=%f, Improvement=%f, Prediction=%f, NA pred=%f\n",
	  dTrainW,
	  dImprovement,
	  dPrediction,
	  (pMissingNode == NULL ? 0.0 : pMissingNode->dPrediction));

  for(i=0; i< cIndent; i++) Rprintf("  ");
  Rprintf("V%d < %f\n",
	  iSplitVar,
	  dSplitValue);
  pLeftNode->PrintSubtree(cIndent+1);
  
  for(i=0; i< cIndent; i++) Rprintf("  ");
  Rprintf("V%d > %f\n",
	  iSplitVar,
	  dSplitValue);
  pRightNode->PrintSubtree(cIndent+1);

  for(i=0; i< cIndent; i++) Rprintf("  ");
  Rprintf("missing\n");
  pMissingNode->PrintSubtree(cIndent+1);
}



signed char CNodeContinuous::WhichNode
(
    CDataset *pData,
    unsigned long iObs
)
{
    signed char ReturnValue = 0;
    double dX = pData->adX[iSplitVar*(pData->cRows) + iObs];

    if(!ISNA(dX))
    {
        if(dX < dSplitValue)
        {
            ReturnValue = -1;
        }
        else
        {
            ReturnValue = 1;
        }
    }
    // if missing value returns 0

    return ReturnValue;
}


signed char CNodeContinuous::WhichNode
(
    double *adX,
    unsigned long cRow,
    unsigned long cCol,
    unsigned long iRow
)
{
    signed char ReturnValue = 0;
    double dX = adX[iSplitVar*cRow + iRow];

    if(!ISNA(dX))
    {
        if(dX < dSplitValue)
        {
            ReturnValue = -1;
        }
        else
        {
            ReturnValue = 1;
        }
    }
    // if missing value returns 0

    return ReturnValue;
}



void CNodeContinuous::RecycleSelf
(
 CNodeFactory *pNodeFactory
)
{
  pNodeFactory->RecycleNode(this);
    
};



void CNodeContinuous::TransferTreeToRList
(
    int &iNodeID,
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
    double dShrinkage
)
{
  int iThisNodeID = iNodeID;
  
  aiSplitVar[iThisNodeID] = iSplitVar;
  adSplitPoint[iThisNodeID] = dSplitValue;
  adErrorReduction[iThisNodeID] = dImprovement;
  adWeight[iThisNodeID] = dTrainW;
  adPred[iThisNodeID] = dShrinkage*dPrediction;
  
  
  iNodeID++;
  aiLeftNode[iThisNodeID] = iNodeID;
  pLeftNode->TransferTreeToRList(iNodeID,
				 pData,
				 aiSplitVar,
				 adSplitPoint,
				 aiLeftNode,
				 aiRightNode,
				 aiMissingNode,
				 adErrorReduction,
				 adWeight,
				 adPred,
				 vecSplitCodes,
				 cCatSplitsOld,
				 dShrinkage);

  aiRightNode[iThisNodeID] = iNodeID;
  pRightNode->TransferTreeToRList(iNodeID,
				  pData,
				  aiSplitVar,
				  adSplitPoint,
				  aiLeftNode,
				  aiRightNode,
				  aiMissingNode,
				  adErrorReduction,
				  adWeight,
				  adPred,
				  vecSplitCodes,
				  cCatSplitsOld,
				  dShrinkage);

  aiMissingNode[iThisNodeID] = iNodeID;
  pMissingNode->TransferTreeToRList(iNodeID,
				    pData,
				    aiSplitVar,
				    adSplitPoint,
				    aiLeftNode,
				    aiRightNode,
				    aiMissingNode,
				    adErrorReduction,
				    adWeight,
				    adPred,
				    vecSplitCodes,
				    cCatSplitsOld,
				    dShrinkage);
}


