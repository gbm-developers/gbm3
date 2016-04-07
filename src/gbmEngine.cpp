//  GBM by Greg Ridgeway  Copyright (C) 2003
//#define NOISY_DEBUG
#include <algorithm>
#include "gbmEngine.h"

CGBM::CGBM(configStructs GBMParams)
{
	// Set up checks for initialization
    fInitialized = false;

    // Create and Initialize dist
    pDataCont = new CGBMDataContainer(GBMParams.GetDataConfig());
    pDataCont->Initialize();

    // Create Tree Components
    pTreeComp = new CTreeComps(GBMParams.GetTreeConfig());
	adZ.assign(pDataCont->getData()->nrow(), 0);
	fInitialized = true;
}


CGBM::~CGBM()
{
	delete pDataCont;
	delete pTreeComp;
}

void CGBM::FitLearner
(
  double *adF,
  double &dTrainError,
  double &dValidError,
  double &dOOBagImprove
)
{
  if(!fInitialized)
  {
    throw GBM::failure("GBM not initialized");
  }

  dTrainError = 0.0;
  dValidError = 0.0;
  dOOBagImprove = 0.0;

  // Initialize adjustments to function estimate
  std::vector<double> adFadj(pDataCont->getData()->nrow(), 0);

  // Bag data
  pDataCont->BagData();

#ifdef NOISY_DEBUG
  Rprintf("Compute working response\n");
#endif

  // Compute Residuals and fit tree
  pDataCont->ComputeResiduals(&adF[0], &adZ[0]);
  pTreeComp->GrowTrees(pDataCont->getData(), &adZ[0], &adFadj[0]);

  // Now I have adF, adZ, and vecpTermNodes (new node assignments)
  // Fit the best constant within each terminal node
#ifdef NOISY_DEBUG
  Rprintf("fit best constant\n");
#endif

  // Adjust terminal node predictions and shrink
  pDataCont->ComputeBestTermNodePreds(&adF[0], &adZ[0], pTreeComp);
  pTreeComp->AdjustAndShrink(&adFadj[0]);

  // Compute the error improvement within bag
  dOOBagImprove = pDataCont->ComputeBagImprovement(&adF[0], pTreeComp->ShrinkageConstant(), &adFadj[0]);

  // Update the function estimate
  unsigned long i = 0;
  for(i=0; i < pDataCont->getData()->get_trainSize(); i++)
  {
    adF[i] += pTreeComp->ShrinkageConstant() * adFadj[i];

  }

  // Make validation predictions
  dTrainError = pDataCont->ComputeDeviance(&adF[0], false);
  pTreeComp->PredictValid(pDataCont->getData(), &adFadj[0]);

  for(i=pDataCont->getData()->get_trainSize();
	  i < pDataCont->getData()->get_trainSize()+pDataCont->getData()->GetValidSize();
	  i++)
  {
    adF[i] += adFadj[i];
  }
  dValidError = pDataCont->ComputeDeviance(&adF[0], true);

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

const long CGBM::SizeOfFittedTree() const
{
	return pTreeComp->GetSizeOfTree();
}

double CGBM::InitF()
{
	return pDataCont->InitialFunctionEstimate();
}
