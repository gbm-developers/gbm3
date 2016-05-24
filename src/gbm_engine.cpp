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
  double* kFuncEstimate,
  double &trainingerror,
  double &validationerror,
  double &outofbag_improvement
)
{

  trainingerror = 0.0;
  validationerror = 0.0;
  outofbag_improvement = 0.0;

  // Initialize adjustments to function estimate
  std::vector<double> delta_estimates(datacontainer_.get_data().nrow(), 0);

  // Bag data
  datacontainer_.BagData();

#ifdef NOISY_DEBUG
  Rprintf("Compute working response\n");
#endif

  // Compute Residuals and fit tree
  datacontainer_.ComputeResiduals(&kFuncEstimate[0], &residuals_[0]);
  treecomponents_.GrowTrees(datacontainer_.get_data(), &residuals_[0], &delta_estimates[0]);



  // Now I have adF, adZ, and vecpTermNodes (new node assignments)
  // Fit the best constant within each terminal node
#ifdef NOISY_DEBUG
  Rprintf("fit best constant\n");
#endif

  // Adjust terminal node predictions and shrink
  datacontainer_.ComputeBestTermNodePreds(&kFuncEstimate[0], &residuals_[0], treecomponents_);
  treecomponents_.AdjustAndShrink(&delta_estimates[0]);

  // Compute the error improvement within bag
  outofbag_improvement = datacontainer_.ComputeBagImprovement(&kFuncEstimate[0],
						 treecomponents_.get_shrinkage_factor(),
						 &delta_estimates[0]);

  // Update the function estimate
  unsigned long i = 0;
  for(i=0; i < datacontainer_.get_data().get_trainsize(); i++)
  {
    kFuncEstimate[i] += treecomponents_.get_shrinkage_factor() * delta_estimates[i];

  }

  // Make validation predictions
  trainingerror = datacontainer_.ComputeDeviance(&kFuncEstimate[0], false);
  treecomponents_.PredictValid(datacontainer_.get_data(), &delta_estimates[0]);

  for(i=datacontainer_.get_data().get_trainsize();
      i < datacontainer_.get_data().get_trainsize()+datacontainer_.get_data().get_validsize();
      i++)
    {
      kFuncEstimate[i] += delta_estimates[i];
    }

  validationerror = datacontainer_.ComputeDeviance(&kFuncEstimate[0], true);

}


void CGBM::GBMTransferTreeToRList
(
 int* splitvar,
 double* splitvalues,
 int* leftnodes,
 int* rightnodes,
 int* missingnodes,
 double* error_reduction,
 double* weights,
 double* predictions,
 VEC_VEC_CATEGORIES &splitcodes_vec,
 int prev_categorical_splits
 )
{
	treecomponents_.TransferTreeToRList(datacontainer_.get_data(),
				     splitvar,
				     splitvalues,
				     leftnodes,
				     rightnodes,
				     missingnodes,
				     error_reduction,
				     weights,
				     predictions,
				     splitcodes_vec,
				     prev_categorical_splits);
}
