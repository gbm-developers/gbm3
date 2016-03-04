//  GBM by Greg Ridgeway  Copyright (C) 2003
//#define NOISY_DEBUG
#include <algorithm>
#include "gbm_engine.h"

CGBM::CGBM(CDistribution* DistPtr, double dLambda,
	    unsigned long cTrain,
	    unsigned long cFeatures,
	    double dBagFraction,
	    unsigned long cDepth,
	    unsigned long cMinObsInNode,
	    int cGroups )
{
    fInitialized = false;
    if(!DistPtr)
    {
       throw GBM::invalid_argument("GBM object could not be initialized - distribution is null");
    }
    pDist=DistPtr;
    pTreeComp = new CTreeComps(dLambda, cTrain, cFeatures,dBagFraction,
    	  cDepth, cMinObsInNode, cGroups);
}


CGBM::~CGBM()
{
	delete pTreeComp;
}


void CGBM::Initialize()
{
  pNodeFactory.reset(new CNodeFactory());
  pNodeFactory->Initialize(pTreeComp->GetDepth());
  pTreeComp -> Initialize(pDist, pNodeFactory.get());
  fInitialized = true;
}

void CGBM::iterate
(
  double *adF,
  double &dTrainError,
  double &dValidError,
  double &dOOBagImprove,
  int &cNodes
)
{
  if(!fInitialized)
  {
    throw GBM::failure("GBM not initialized");
  }

  dTrainError = 0.0;
  dValidError = 0.0;
  dOOBagImprove = 0.0;

  pTreeComp->AssignTermNodes();
  pTreeComp->BagData(IsPairwise(), pDist);

#ifdef NOISY_DEBUG
  Rprintf("Compute working response\n");
#endif

  pDist->ComputeWorkingResponse(adF,
                               	pTreeComp->GetGrad(),
                                pTreeComp->GetBag(),
                                pTreeComp->GetTrainNo());

  pTreeComp->GrowTrees(pDist, cNodes);

  // Now I have adF, adZ, and vecpTermNodes (new node assignments)
  // Fit the best constant within each terminal node
#ifdef NOISY_DEBUG
  Rprintf("fit best constant\n");
#endif

  pDist->FitBestConstant(&adF[0],
                         pTreeComp->GetGrad(),
                         pTreeComp->GetNodeAssign(),
                         pTreeComp->GetTrainNo(),
                         pTreeComp->GetTermNodes(),
                         (2*cNodes+1)/3, // number of terminal nodes
                         pTreeComp->GetMinNodeObs(),
                         pTreeComp->GetBag(),
                         pTreeComp->GetRespAdj());

  pTreeComp->AdjustAndShrink();
  // update training predictions
  // fill in missing nodes where N < cMinObsInNode


  dOOBagImprove = pDist->BagImprovement(&adF[0],
                                        pTreeComp->GetRespAdj(),
                                        pTreeComp->GetBag(),
                                        pTreeComp->GetLambda(),
                                        pTreeComp->GetTrainNo());
    
  // update the training predictions
  unsigned long i = 0;
  for(i=0; i < pTreeComp->GetTrainNo(); i++)
  {
    adF[i] += pTreeComp->GetLambda() * pTreeComp->RespAdjElem(i);
  }

  dTrainError = pDist->Deviance(adF, pTreeComp->GetTrainNo());

  // update the validation predictions
  pTreeComp->PredictValid(pDist);

  for(i=pTreeComp->GetTrainNo(); i < pTreeComp->GetTrainNo()+pTreeComp->GetValidNo(); i++)
  {
    adF[i] += pTreeComp->RespAdjElem(i);
  }

  dValidError = pDist->Deviance(adF + pTreeComp->GetTrainNo(), pTreeComp->GetValidNo(), true);

}


void CGBM::TransferTreeToRList
(
 int *aiSplitVar,
 double *adSplitPoint,
 int *aiLeftNode,
 int *aiRightNode,
 int *aiMissingNode,
 double *adErrorReduction,
 double *adWeight,
 double *adPred,
 VEC_VEC_CATEGORIES &vecSplitCodes,
 int cCatSplitsOld
 )
{
	pTreeComp->TransferTreeToRList(*(pDist->data_ptr()),
				 aiSplitVar,
				 adSplitPoint,
				 aiLeftNode,
				 aiRightNode,
				 aiMissingNode,
				 adErrorReduction,
				 adWeight,
				 adPred,
				 vecSplitCodes,
				 cCatSplitsOld);
}


