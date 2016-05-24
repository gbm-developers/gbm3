//  GBM by Greg Ridgeway  Copyright (C) 2003
//#define NOISY_DEBUG
#include <algorithm>
#include "gbm_engine.h"

CGBM::CGBM(ConfigStructs& GBMParams) :
  datacontainer_(GBMParams.get_data_config()),
  treecomponents_(GBMParams.get_tree_config()),
  residuals_(datacontainer_.get_data().nrow(), 0) {}


CGBM::~CGBM()
{
}

void CGBM::FitLearner
(
  double *adF,
  double &trainingerror,
  double &validationerror,
  double &outofbag_improvement
)
{

  trainingerror = 0.0;
  validationerror = 0.0;
  outofbag_improvement = 0.0;

  // Initialize adjustments to function estimate
  std::vector<double> adFadj(datacontainer_.get_data().nrow(), 0);

  // Bag data
  datacontainer_.BagData();

#ifdef NOISY_DEBUG
  Rprintf("Compute working response\n");
#endif

  // Compute Residuals and fit tree
  datacontainer_.ComputeResiduals(&adF[0], &residuals_[0]);
  treecomponents_.GrowTrees(datacontainer_.get_data(), &residuals_[0], &adFadj[0]);



  // Now I have adF, adZ, and vecpTermNodes (new node assignments)
  // Fit the best constant within each terminal node
#ifdef NOISY_DEBUG
  Rprintf("fit best constant\n");
#endif

  // Adjust terminal node predictions and shrink
  datacontainer_.ComputeBestTermNodePreds(&adF[0], &residuals_[0], treecomponents_);
  treecomponents_.AdjustAndShrink(&adFadj[0]);

  // Compute the error improvement within bag
  outofbag_improvement = datacontainer_.ComputeBagImprovement(&adF[0],
						 treecomponents_.get_shrinkage_factor(),
						 &adFadj[0]);

  // Update the function estimate
  unsigned long i = 0;
  for(i=0; i < datacontainer_.get_data().get_trainsize(); i++)
  {
    adF[i] += treecomponents_.get_shrinkage_factor() * adFadj[i];

  }

  // Make validation predictions
  trainingerror = datacontainer_.ComputeDeviance(&adF[0], false);
  treecomponents_.PredictValid(datacontainer_.get_data(), &adFadj[0]);

  for(i=datacontainer_.get_data().get_trainsize();
      i < datacontainer_.get_data().get_trainsize()+datacontainer_.get_data().get_validsize();
      i++)
    {
      adF[i] += adFadj[i];
    }

  validationerror = datacontainer_.ComputeDeviance(&adF[0], true);

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
	treecomponents_.TransferTreeToRList(datacontainer_.get_data(),
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
