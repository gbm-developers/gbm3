//  GBM by Greg Ridgeway  Copyright (C) 2003
//#define NOISY_DEBUG
#include <algorithm>
#include "gbmEngine.h"

CGBM::CGBM()
{
	// Set up checks for initialization
    fInitialized = false;
    hasDataAndDist = false;
    hasTreeContainer = false;

    // Containers
    pDataCont = NULL;
    pTreeComp = NULL;
    pNodeFactory = new CNodeFactory();
}


CGBM::~CGBM()
{
	delete pDataCont;
	delete pTreeComp;
	delete pNodeFactory;
}

void CGBM::SetDataAndDistribution(SEXP radY, SEXP radOffset, SEXP radX, SEXP raiXOrder,
        SEXP radWeight, SEXP racVarClasses,
        SEXP ralMonotoneVar, SEXP radMisc, const std::string& family,
		const int cTrain, int& cGroups)
{

	pDataCont=new CGBMDataContainer(radY, radOffset, radX, raiXOrder,
        radWeight, racVarClasses, ralMonotoneVar, radMisc, family, cTrain, cGroups);
	hasDataAndDist = true;
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
	hasTreeContainer = true;
}

void CGBM::Initialize()
{
	// Throw error if initialization called incorrectly
	if(!(hasDataAndDist && hasTreeContainer))
	{
		throw GBM::failure("GBM object could not be built - missing: "
				"data or distribution or weak learner");
	}
	pNodeFactory->NodeFactoryInitialize(pTreeComp->GetDepth());
	pDataCont-> Initialize();
	pTreeComp -> TreeInitialize(pDataCont->getData(), pNodeFactory);
	fInitialized = true;
}

void CGBM::Iterate
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
  pTreeComp->BagData(IsPairwise(), pDataCont->getDist());

#ifdef NOISY_DEBUG
  Rprintf("Compute working response\n");
#endif

  pDataCont->ComputeResiduals(&adF[0], pTreeComp);
  pTreeComp->GrowTrees(pDataCont->getData(), cNodes);

  // Now I have adF, adZ, and vecpTermNodes (new node assignments)
  // Fit the best constant within each terminal node
#ifdef NOISY_DEBUG
  Rprintf("fit best constant\n");
#endif

  pDataCont->ComputeBestTermNodePreds(&adF[0], pTreeComp, cNodes);
  pTreeComp->AdjustAndShrink();

  // update training predictions
  // fill in missing nodes where N < cMinObsInNode
  dOOBagImprove = pDataCont->ComputeBagImprovement(&adF[0], pTreeComp);

  // update the training predictions
  unsigned long i = 0;
  for(i=0; i < pTreeComp->GetTrainNo(); i++)
  {
    adF[i] += pTreeComp->GetLambda() * pTreeComp->RespAdjElem(i);
  }
  dTrainError = pDataCont->ComputeDeviance(&adF[0], pTreeComp, false);

  // update the validation predictions
  pTreeComp->PredictValid(pDataCont->getData());
  for(i=pTreeComp->GetTrainNo(); i < pTreeComp->GetTrainNo()+pTreeComp->GetValidNo(); i++)
  {
    adF[i] += pTreeComp->RespAdjElem(i);
  }

  dValidError = pDataCont->ComputeDeviance(&adF[0], pTreeComp, true);

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
	pTreeComp->TransferTreeToRList(*(pDataCont->getData()),
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
	pDataCont->InitializeFunctionEstimate(dInitF, cLength);
}

void CGBM::UpdateParams(const double *adF,
        			      unsigned long cLength)
{
	pDataCont->UpdateParams(adF, cLength);
}
