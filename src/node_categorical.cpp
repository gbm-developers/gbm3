//  GBM by Greg Ridgeway  Copyright (C) 2003

#include "node_categorical.h"
#include "node_factory.h"


CNodeCategorical::~CNodeCategorical()
{
    #ifdef NOISY_DEBUG
    Rprintf("categorical destructor\n");
    #endif
}


void CNodeCategorical::PrintSubtree
(
    unsigned long cIndent
)
{
  unsigned long i = 0;
  const std::size_t cLeftCategory = aiLeftCategory.size();
  
  for(i=0; i< cIndent; i++) Rprintf("  ");
  Rprintf("N=%f, Improvement=%f, Prediction=%f, NA pred=%f\n",
	  dTrainW,
	  dImprovement,
	  dPrediction,
	  (pMissingNode == NULL ? 0.0 : pMissingNode->dPrediction));

  for(i=0; i< cIndent; i++) Rprintf("  ");
  Rprintf("V%d in ",iSplitVar);
  for(i=0; i<cLeftCategory; i++)
    {
      Rprintf("%d",aiLeftCategory[i]);
      if(i<cLeftCategory-1) Rprintf(",");
    }
  Rprintf("\n");
  pLeftNode->PrintSubtree(cIndent+1);

  for(i=0; i< cIndent; i++) Rprintf("  ");
  Rprintf("V%d not in ",iSplitVar);
  for(i=0; i<cLeftCategory; i++)
    {
      Rprintf("%d",aiLeftCategory[i]);
      if(i<cLeftCategory-1) Rprintf(",");
    }
  Rprintf("\n");
  pRightNode->PrintSubtree(cIndent+1);
  
  for(i=0; i< cIndent; i++) Rprintf("  ");
  Rprintf("missing\n");
  pMissingNode->PrintSubtree(cIndent+1);
}


signed char CNodeCategorical::WhichNode
(
    CDataset &data,
    unsigned long iObs
)
{
    signed char ReturnValue = 0;
    double dX = data.x_value(iObs, iSplitVar);

    if(!ISNA(dX))
    {
      if(std::find(aiLeftCategory.begin(),
		   aiLeftCategory.end(),
		   (ULONG)dX) != aiLeftCategory.end())
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


signed char CNodeCategorical::WhichNode
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
      if(std::find(aiLeftCategory.begin(),
		   aiLeftCategory.end(),
		   (ULONG)dX) != aiLeftCategory.end())
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




void CNodeCategorical::RecycleSelf(CNodeFactory *pNodeFactory) {
  pNodeFactory->RecycleNode(this);
}



void CNodeCategorical::TransferTreeToRList
(
 int &iNodeID,
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
 double dShrinkage
)
{
  
  int iThisNodeID = iNodeID;
  unsigned long cCatSplits = vecSplitCodes.size();
  unsigned long i = 0;
  int cLevels = data.varclass_ptr()[iSplitVar];
  const std::size_t cLeftCategory = aiLeftCategory.size();
  
  aiSplitVar[iThisNodeID] = iSplitVar;
  adSplitPoint[iThisNodeID] = cCatSplits+cCatSplitsOld; // 0 based
  adErrorReduction[iThisNodeID] = dImprovement;
  adWeight[iThisNodeID] = dTrainW;
  adPred[iThisNodeID] = dShrinkage*dPrediction;
  
  vecSplitCodes.push_back(VEC_CATEGORIES());
  
  vecSplitCodes[cCatSplits].resize(cLevels,1);
  for(i=0; i<cLeftCategory; i++)
    {
      vecSplitCodes[cCatSplits][aiLeftCategory[i]] = -1;
    }

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


