//  GBM by Greg Ridgeway  Copyright (C) 2003
//#define NOISY_DEBUG
#include <algorithm>
#include "gbm_engine.h"

CGBMEngine::CGBMEngine(DataDistParams& datadistparams,
		TreeParams& treeparams)
    : datacontainer_(datadistparams),
      tree_params_(treeparams),
      residuals_(datacontainer_.get_data().nrow(), 0) {}

CGBMEngine::~CGBMEngine() {}

FittedLearner* CGBMEngine::FitLearner(double* func_estimate) {

  // Initialize adjustments to function estimate
  std::vector<double> delta_estimates(datacontainer_.get_data().nrow(), 0);

  // Bag data
  datacontainer_.BagData();

  // Set up tree
  std::auto_ptr<CCARTTree> tree(new CCARTTree(tree_params_));

  // Compute Residuals and fit tree
  datacontainer_.ComputeResiduals(&func_estimate[0], residuals_);
  tree->Grow(residuals_, datacontainer_.get_data(),
		  datacontainer_.get_bag(), delta_estimates);

  // Now I have adF, adZ, and vecpTermNodes (new node assignments)
  // Fit the best constant within each terminal node

  // Adjust terminal node predictions and shrink
  datacontainer_.ComputeBestTermNodePreds(&func_estimate[0], residuals_,
                                          *tree.get());
  tree->Adjust(delta_estimates);

  // Compute the error improvement within bag
  double oobag_improv = datacontainer_.ComputeBagImprovement(
      &func_estimate[0], tree->get_shrinkage_factor(), delta_estimates);

  // Update the function estimate
  unsigned long i = 0;
  for (i = 0; i < datacontainer_.get_data().get_trainsize(); i++) {
    func_estimate[i] += tree->get_shrinkage_factor() * delta_estimates[i];
  }

  // Make validation predictions
  double train_error = datacontainer_.ComputeDeviance(&func_estimate[0], false);
  tree->PredictValid(datacontainer_.get_data(),
                     datacontainer_.get_data().get_validsize(),
                     delta_estimates);

  for (i = datacontainer_.get_data().get_trainsize();
       i < datacontainer_.get_data().get_trainsize() +
               datacontainer_.get_data().get_validsize();
       i++) {
    func_estimate[i] += delta_estimates[i];
  }

  double valid_error = datacontainer_.ComputeDeviance(&func_estimate[0], true);
  std::auto_ptr<FittedLearner> fit(new FittedLearner(tree, datacontainer_.get_data(),
		  	  	  	  train_error, valid_error, oobag_improv));
  return fit.release();
}
