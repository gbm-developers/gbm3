//  GBM by Greg Ridgeway  Copyright (C) 2003
//#define NOISY_DEBUG
#include <algorithm>
#include "gbm_engine.h"

CGBM::CGBM()
{
    fInitialized = false;
    pDist = NULL;
    pTreeComp = NULL;
    pNodeFactory = new CNodeFactory();
}


CGBM::~CGBM()
{
	delete pDist;
	delete pTreeComp;
	delete pNodeFactory;
}

void CGBM::SetDataAndDistribution(const CDataset& data, SEXP radMisc, const std::string& family,
		const int cTrain, int& cGroups)
{
	pDist=gbm_setup(data, radMisc, family, cTrain, cGroups);
}

void CGBM::SetTreeContainer(double dLambda,
   	    unsigned long cTrain,
   	    unsigned long cFeatures,
   	    double dBagFraction,
   	    unsigned long cDepth,
   	    unsigned long cMinObsInNode,
   	    int cGroups)
{
	pTreeComp = new CTreeComps(dLambda, cTrain, cFeatures,dBagFraction,
	    	  cDepth, cMinObsInNode, cGroups);
}

void CGBM::Initialize()
{
  pNodeFactory->NodeFactoryInitialize(pTreeComp->GetDepth());
  pDist-> Initialize();
  pTreeComp -> TreeInitialize(pDist, pNodeFactory);
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


void CGBM::GBMTransferTreeToRList
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

void CGBM::InitF(double &dInitF, unsigned long cLength)
{
	pDist->InitF(dInitF, cLength);
}

void CGBM::UpdateParams(const double *adF,
        			      unsigned long cLength)
{
	pDist->UpdateParams(adF, cLength);
}
