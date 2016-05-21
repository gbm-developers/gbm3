//  GBM by Greg Ridgeway  Copyright (C) 2003
//#define NOISY_DEBUG
#include <algorithm>
#include "gbmEngine.h"

CGBM::CGBM(configStructs& GBMParams) :
  dataCont(GBMParams.GetDataConfig()),
  treeComp(GBMParams.GetTreeConfig()),
  adZ(dataCont.getData().nrow(), 0) {}


CGBM::~CGBM()
{
}

void CGBM::FitLearner
(
  double *adF,
  double &dTrainError,
  double &dValidError,
  double &dOOBagImprove
)
{

  dTrainError = 0.0;
  dValidError = 0.0;
  dOOBagImprove = 0.0;

  // Initialize adjustments to function estimate
  std::vector<double> adFadj(dataCont.getData().nrow(), 0);

  // Bag data
  dataCont.BagData();

#ifdef NOISY_DEBUG
  Rprintf("Compute working response\n");
#endif

  // Compute Residuals and fit tree
  dataCont.ComputeResiduals(&adF[0], &adZ[0]);
  treeComp.GrowTrees(dataCont.getData(), &adZ[0], &adFadj[0]);



  // Now I have adF, adZ, and vecpTermNodes (new node assignments)
  // Fit the best constant within each terminal node
#ifdef NOISY_DEBUG
  Rprintf("fit best constant\n");
#endif

  // Adjust terminal node predictions and shrink
  dataCont.ComputeBestTermNodePreds(&adF[0], &adZ[0], treeComp);
  treeComp.AdjustAndShrink(&adFadj[0]);

  // Compute the error improvement within bag
  dOOBagImprove = dataCont.ComputeBagImprovement(&adF[0],
						 treeComp.ShrinkageConstant(),
						 &adFadj[0]);

  // Update the function estimate
  unsigned long i = 0;
  for(i=0; i < dataCont.getData().get_trainSize(); i++)
  {
    adF[i] += treeComp.ShrinkageConstant() * adFadj[i];

  }

  // Make validation predictions
  dTrainError = dataCont.ComputeDeviance(&adF[0], false);
  treeComp.PredictValid(dataCont.getData(), &adFadj[0]);

  for(i=dataCont.getData().get_trainSize();
      i < dataCont.getData().get_trainSize()+dataCont.getData().GetValidSize();
      i++)
    {
      adF[i] += adFadj[i];
    }

  dValidError = dataCont.ComputeDeviance(&adF[0], true);

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
	treeComp.TransferTreeToRList(dataCont.getData(),
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
  return treeComp.GetSizeOfTree();
}

double CGBM::InitF()
{
  return dataCont.InitialFunctionEstimate();
}
