//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//
//  File:       node_search.cpp
//
//------------------------------------------------------------------------------
#include <cassert>
#include "node_search.h"

CNodeSearch::CNodeSearch()
{
    iBestSplitVar = 0;

    dBestSplitValue = 0.0;
    fIsSplit = false;

    dBestMissingTotalW = 0.0;
    dCurrentMissingTotalW = 0.0;
    dBestMissingSumZ = 0.0;
    dCurrentMissingSumZ = 0.0;

    adGroupSumZ.resize(1024);
    adGroupW.resize(1024);
    acGroupN.resize(1024);
    adGroupMean.resize(1024);
    aiCurrentCategory.resize(1024);
    aiBestCategory.resize(1024);

    iRank = UINT_MAX;
}


CNodeSearch::~CNodeSearch()
{
}


void CNodeSearch::Initialize
(
    unsigned long cMinObsInNode
)
{
    this->cMinObsInNode = cMinObsInNode;
}


void CNodeSearch::IncorporateObs
(
    double dX,
    double dZ,
    double dW,
    long lMonotone
)
{
    static double dWZ = 0.0;

    if(fIsSplit) return;

    dWZ = dW*dZ;

    if(ISNA(dX))
    {
        dCurrentMissingSumZ += dWZ;
        dCurrentMissingTotalW += dW;
        cCurrentMissingN++;
        dCurrentRightSumZ -= dWZ;
        dCurrentRightTotalW -= dW;
        cCurrentRightN--;
    }
    else if(cCurrentVarClasses == 0)   // variable is continuous
    {
        if(dLastXValue > dX)
        {
	  throw GBM::failure("Observations are not in order. gbm() was unable to build an index for the design matrix. Could be a bug in gbm or an unusual data type in data.");
	}

        // Evaluate the current split
        // the newest observation is still in the right child
        dCurrentSplitValue = 0.5*(dLastXValue + dX);
        if((dLastXValue != dX) &&
            (cCurrentLeftN >= cMinObsInNode) &&
            (cCurrentRightN >= cMinObsInNode) &&
            ((lMonotone==0) ||
            (lMonotone*(dCurrentRightSumZ*dCurrentLeftTotalW -
                        dCurrentLeftSumZ*dCurrentRightTotalW) > 0)))
        {
            dCurrentImprovement =
                CNode::Improvement(dCurrentLeftTotalW,dCurrentRightTotalW,
                                    dCurrentMissingTotalW,
                                    dCurrentLeftSumZ,dCurrentRightSumZ,
                                    dCurrentMissingSumZ);
            if(dCurrentImprovement > dBestImprovement)
            {
                iBestSplitVar = iCurrentSplitVar;
                dBestSplitValue = dCurrentSplitValue;
                cBestVarClasses = 0;

                dBestLeftSumZ    = dCurrentLeftSumZ;
                dBestLeftTotalW  = dCurrentLeftTotalW;
                cBestLeftN       = cCurrentLeftN;
                dBestRightSumZ   = dCurrentRightSumZ;
                dBestRightTotalW = dCurrentRightTotalW;
                cBestRightN      = cCurrentRightN;
                dBestImprovement = dCurrentImprovement;
            }
        }

        // now move the new observation to the left
        // if another observation arrives we will evaluate this
        dCurrentLeftSumZ += dWZ;
        dCurrentLeftTotalW += dW;
        cCurrentLeftN++;
        dCurrentRightSumZ -= dWZ;
        dCurrentRightTotalW -= dW;
        cCurrentRightN--;

        dLastXValue = dX;
    }
    else // variable is categorical, evaluates later
    {
        adGroupSumZ[(unsigned long)dX] += dWZ;
        adGroupW[(unsigned long)dX] += dW;
        acGroupN[(unsigned long)dX] ++;
    }
}



void CNodeSearch::Set
(
    double dSumZ,
    double dTotalW,
    unsigned long cTotalN,
    CNodeTerminal *pThisNode,
    CNode **ppParentPointerToThisNode,
    CNodeFactory *pNodeFactory
)
{
    dInitSumZ = dSumZ;
    dInitTotalW = dTotalW;
    cInitN = cTotalN;

    dBestLeftSumZ       = 0.0;
    dBestLeftTotalW     = 0.0;
    cBestLeftN          = 0;
    dCurrentLeftSumZ    = 0.0;
    dCurrentLeftTotalW  = 0.0;
    cCurrentLeftN       = 0;

    dBestRightSumZ      = dSumZ;
    dBestRightTotalW    = dTotalW;
    cBestRightN         = cTotalN;
    dCurrentRightSumZ   = 0.0;
    dCurrentRightTotalW = dTotalW;
    cCurrentRightN      = cTotalN;

    dBestMissingSumZ      = 0.0;
    dBestMissingTotalW    = 0.0;
    cBestMissingN         = 0;
    dCurrentMissingSumZ   = 0.0;
    dCurrentMissingTotalW = 0.0;
    cCurrentMissingN      = 0;

    dBestImprovement    = 0.0;
    iBestSplitVar       = UINT_MAX;

    dCurrentImprovement = 0.0;
    iCurrentSplitVar    = UINT_MAX;
    dCurrentSplitValue  = -HUGE_VAL;

    fIsSplit = false;

    this->pThisNode = pThisNode;
    this->ppParentPointerToThisNode = ppParentPointerToThisNode;
    this->pNodeFactory = pNodeFactory;
}


void CNodeSearch::ResetForNewVar
(
    unsigned long iWhichVar,
    long cCurrentVarClasses
)
{
  if(fIsSplit) return;

  assert(cCurrentVarClasses <= adGroupSumZ.size());

  std::fill(adGroupSumZ.begin(), adGroupSumZ.begin() + cCurrentVarClasses, 0);
  std::fill(adGroupW.begin(), adGroupW.begin() + cCurrentVarClasses, 0);
  std::fill(acGroupN.begin(), acGroupN.begin() + cCurrentVarClasses, 0);
  
  iCurrentSplitVar = iWhichVar;
  this->cCurrentVarClasses = cCurrentVarClasses;

  dCurrentLeftSumZ      = 0.0;
  dCurrentLeftTotalW    = 0.0;
  cCurrentLeftN         = 0;
  dCurrentRightSumZ     = dInitSumZ;
  dCurrentRightTotalW   = dInitTotalW;
  cCurrentRightN        = cInitN;
  dCurrentMissingSumZ   = 0.0;
  dCurrentMissingTotalW = 0.0;
  cCurrentMissingN      = 0;
  
  dCurrentImprovement = 0.0;
  
  dLastXValue = -HUGE_VAL;
}



void CNodeSearch::WrapUpCurrentVariable()
{
  if(iCurrentSplitVar == iBestSplitVar)
    {
      if(cCurrentMissingN > 0)
        {
	  dBestMissingSumZ   = dCurrentMissingSumZ;
	  dBestMissingTotalW = dCurrentMissingTotalW;
	  cBestMissingN      = cCurrentMissingN;
        }
      else // DEBUG: consider a weighted average with parent node?
        {
	  dBestMissingSumZ   = dInitSumZ;
	  dBestMissingTotalW = dInitTotalW;
	  cBestMissingN      = 0;
        }
    }
}



void CNodeSearch::EvaluateCategoricalSplit()
{
  long i=0;
  long j=0;
  unsigned long cFiniteMeans = 0;
  
  if(fIsSplit) return;
  
  if(cCurrentVarClasses == 0)
    {
      throw GBM::invalid_argument();
    }

  cFiniteMeans = 0;
  for(i=0; i<cCurrentVarClasses; i++)
    {
      aiCurrentCategory[i] = i;
      if(adGroupW[i] != 0.0)
        {
	  adGroupMean[i] = adGroupSumZ[i]/adGroupW[i];
	  cFiniteMeans++;
        }
      else
        {
	  adGroupMean[i] = HUGE_VAL;
        }
    }
  
  rsort_with_index(&adGroupMean[0],&aiCurrentCategory[0],cCurrentVarClasses);
    
  // if only one group has a finite mean it will not consider
  // might be all are missing so no categories enter here
  for(i=0; (cFiniteMeans>1) && ((ULONG)i<cFiniteMeans-1); i++)
    {
      dCurrentSplitValue = (double)i;
      
      dCurrentLeftSumZ    += adGroupSumZ[aiCurrentCategory[i]];
      dCurrentLeftTotalW  += adGroupW[aiCurrentCategory[i]];
      cCurrentLeftN       += acGroupN[aiCurrentCategory[i]];
      dCurrentRightSumZ   -= adGroupSumZ[aiCurrentCategory[i]];
      dCurrentRightTotalW -= adGroupW[aiCurrentCategory[i]];
      cCurrentRightN      -= acGroupN[aiCurrentCategory[i]];
      
      dCurrentImprovement =
	CNode::Improvement(dCurrentLeftTotalW,dCurrentRightTotalW,
			   dCurrentMissingTotalW,
			   dCurrentLeftSumZ,dCurrentRightSumZ,
			   dCurrentMissingSumZ);
      if((cCurrentLeftN >= cMinObsInNode) &&
	 (cCurrentRightN >= cMinObsInNode) &&
	 (dCurrentImprovement > dBestImprovement))
        {
	  dBestSplitValue = dCurrentSplitValue;
	  if(iBestSplitVar != iCurrentSplitVar)
            {
	      iBestSplitVar = iCurrentSplitVar;
	      cBestVarClasses = cCurrentVarClasses;
	      std::copy(aiCurrentCategory.begin(),
			aiCurrentCategory.end(),
			aiBestCategory.begin());
            }
	  
	  dBestLeftSumZ      = dCurrentLeftSumZ;
	  dBestLeftTotalW    = dCurrentLeftTotalW;
	  cBestLeftN         = cCurrentLeftN;
	  dBestRightSumZ     = dCurrentRightSumZ;
	  dBestRightTotalW   = dCurrentRightTotalW;
	  cBestRightN        = cCurrentRightN;
	  dBestImprovement   = dCurrentImprovement;
        }
    }
}




void CNodeSearch::SetupNewNodes
(
    PCNodeNonterminal &pNewSplitNode,
    PCNodeTerminal &pNewLeftNode,
    PCNodeTerminal &pNewRightNode,
    PCNodeTerminal &pNewMissingNode
)
{
    CNodeContinuous *pNewNodeContinuous = NULL;
    CNodeCategorical *pNewNodeCategorical = NULL;
    unsigned long i=0;

    pNewLeftNode    = pNodeFactory->GetNewNodeTerminal();
    pNewRightNode   = pNodeFactory->GetNewNodeTerminal();
    pNewMissingNode = pNodeFactory->GetNewNodeTerminal();

    // set up a continuous split
    if(cBestVarClasses==0)
    {
        pNewNodeContinuous = pNodeFactory->GetNewNodeContinuous();

        pNewNodeContinuous->dSplitValue = dBestSplitValue;
        pNewNodeContinuous->iSplitVar = iBestSplitVar;

        pNewSplitNode = pNewNodeContinuous;
    }
    else
    {
        // get a new categorical node and its branches
        pNewNodeCategorical = pNodeFactory->GetNewNodeCategorical();

        // set up the categorical split
        pNewNodeCategorical->iSplitVar = iBestSplitVar;
        pNewNodeCategorical->aiLeftCategory.resize(1 + (ULONG)dBestSplitValue);
	std::copy(aiBestCategory.begin(),
		  aiBestCategory.begin() + pNewNodeCategorical->aiLeftCategory.size(),
		  pNewNodeCategorical->aiLeftCategory.begin());

        pNewSplitNode = pNewNodeCategorical;
    }

    *ppParentPointerToThisNode = pNewSplitNode;

    pNewSplitNode->dPrediction  = pThisNode->dPrediction;
    pNewSplitNode->dImprovement = dBestImprovement;
    pNewSplitNode->dTrainW      = pThisNode->dTrainW;
    pNewSplitNode->pLeftNode    = pNewLeftNode;
    pNewSplitNode->pRightNode   = pNewRightNode;
    pNewSplitNode->pMissingNode = pNewMissingNode;

    pNewLeftNode->dPrediction    = dBestLeftSumZ/dBestLeftTotalW;
    pNewLeftNode->dTrainW        = dBestLeftTotalW;
    pNewLeftNode->cN             = cBestLeftN;
    pNewRightNode->dPrediction   = dBestRightSumZ/dBestRightTotalW;
    pNewRightNode->dTrainW       = dBestRightTotalW;
    pNewRightNode->cN            = cBestRightN;
    pNewMissingNode->dPrediction = dBestMissingSumZ/dBestMissingTotalW;
    pNewMissingNode->dTrainW     = dBestMissingTotalW;
    pNewMissingNode->cN          = cBestMissingN;

    pThisNode->RecycleSelf(pNodeFactory);
}
