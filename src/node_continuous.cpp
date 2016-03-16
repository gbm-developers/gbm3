//  GBM by Greg Ridgeway  Copyright (C) 2003

#include "node_continuous.h"


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
    const CDataset &data,
    unsigned long iObs
)
{
    signed char ReturnValue = 0;
    double dX = data.x_value(iObs, iSplitVar);

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

void CNodeContinuous::TransferTreeToRList
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
	if(isTerminal)
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
	else
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
					 data,
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
					  data,
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
						data,
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
  
}


