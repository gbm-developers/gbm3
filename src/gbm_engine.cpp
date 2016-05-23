//  GBM by Greg Ridgeway  Copyright (C) 2003
//#define NOISY_DEBUG
#include <algorithm>
#include "gbm_engine.h"

CGBM::CGBM(ConfigStructs& GBMParams) :
  dataCont(GBMParams.get_data_config()),
  treeComp(GBMParams.get_tree_config()),
  adZ(dataCont.get_data().nrow(), 0) {}


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
  std::vector<double> adFadj(dataCont.get_data().nrow(), 0);

  // Bag data
  dataCont.BagData();

#ifdef NOISY_DEBUG
  Rprintf("Compute working response\n");
#endif

  // Compute Residuals and fit tree
  dataCont.ComputeResiduals(&adF[0], &adZ[0]);
  treeComp.GrowTrees(dataCont.get_data(), &adZ[0], &adFadj[0]);



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
						 treeComp.get_shrinkage_factor(),
						 &adFadj[0]);

  // Update the function estimate
  unsigned long i = 0;
  for(i=0; i < dataCont.get_data().get_trainsize(); i++)
  {
    adF[i] += treeComp.get_shrinkage_factor() * adFadj[i];

  }

  // Make validation predictions
  dTrainError = dataCont.ComputeDeviance(&adF[0], false);
  treeComp.PredictValid(dataCont.get_data(), &adFadj[0]);

  for(i=dataCont.get_data().get_trainsize();
      i < dataCont.get_data().get_trainsize()+dataCont.get_data().get_validsize();
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
	treeComp.TransferTreeToRList(dataCont.get_data(),
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
