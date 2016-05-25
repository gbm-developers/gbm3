//-----------------------------------
//
// File: gbmDataContainer.cpp
//
// Description: class that contains the data and distribution.
//
//-----------------------------------

//------------------------------
// Includes
//------------------------------
#include "gbm_datacontainer.h"
#include "pairwise.h"

//----------------------------------------
// Function Members - Public
//----------------------------------------
//-----------------------------------
// Function: CGBMDataContainer
//
// Returns: none
//
// Description: Default constructor for gbm data container.
//
// Parameters: ...
//-----------------------------------
CGBMDataContainer::CGBMDataContainer(DataDistParams& datadist_config):
        data_(datadist_config)
{
	//Initialize the factory and then use to get the disribution
	distfactory_ = new DistributionFactory();
	distptr_ = distfactory_ -> CreateDist(datadist_config);
	distptr_->Initialize(data_);
}

//-----------------------------------
// Function: ~CGBMDataContainer
//
// Returns: none
//
// Description: Default destructor for GBM Data Container
//
// Parameters: none
//-----------------------------------
CGBMDataContainer::~CGBMDataContainer()
{
	delete distptr_;
	delete distfactory_;
}

//-----------------------------------
// Function: InitializeFunctionEstimate
//
// Returns: none
//
// Description: Initialize the function fit.
//
// Parameters: double& - reference to the initial function estimate (a constant)
//    unsigned long - the number of predictors the fit must provide response estimates for
//
//-----------------------------------
double CGBMDataContainer::InitialFunctionEstimate()
{
  return get_dist()->InitF(data_);
}

//-----------------------------------
// Function: ComputeResiduals
//
// Returns: nonei = 0; i < data.getNumUniquePatient()
//
// Description: Compute the residuals associated with the distributions loss function.
//
// Parameters: const double ptr - ptr to the function estimates for each predictor
//    CTreeComps ptr - ptr to the tree components container in the gbm.
//
//-----------------------------------
void CGBMDataContainer::ComputeResiduals(const double* kFuncEstimate, double* residuals)
{
	get_dist()->ComputeWorkingResponse(data_, kFuncEstimate, residuals);
}

//-----------------------------------
// Function: ComputeBestTermNodePreds
//
// Returns: none
//
// Description: Fit the best constants (predictions) to the terminal nodes.
//
// Parameters: const double ptr - ptr to function estimates for each predictor
//    CTreeComps ptr - ptr to the tree components container in the gbm
//    int& - reference to the number of nodes in the tree.
//-----------------------------------
void CGBMDataContainer::ComputeBestTermNodePreds(const double* kFuncEstimate, double* residuals, CCARTTree& tree)
{
  get_dist()->FitBestConstant(get_data(),
			     &kFuncEstimate[0],
			     (2*tree.size_of_tree()+1)/3, // number of terminal nodes
			     &residuals[0],
			     tree);
}

//-----------------------------------
// Function: ComputeDeviance
//
// Returns: double
//
// Description: Compute the deviance (error) of the fit on training/validation data.
//
// Parameters: const double ptr - ptr to function estimates for each predictor
//    CTreeComps ptr - ptr to the tree components container in the gbm
//    bool - bool which indicates whether it is the training or validation data used.
//
//-----------------------------------
double CGBMDataContainer::ComputeDeviance(const double* kFuncEstimate, bool is_validationset)
{
  if(!(is_validationset))
    {
      return get_dist()->Deviance(data_, kFuncEstimate);
    }
  else
    {
      return get_dist()->Deviance(data_, kFuncEstimate + data_.get_trainsize(), true);
    }
}

//-----------------------------------
// Function: ComputeBagImprovement
//
// Returns: double
//
// Description: Compute the improvement from combining decision trees
//
// Parameters: const double ptr - ptr to function estimates for each predictor
//    CTreeComps ptr - ptr to the tree components container in the gbm
//
//-----------------------------------
double CGBMDataContainer::ComputeBagImprovement(const double* kFuncEstimate, const double kShrinkage, const double* kDeltaEstimate)
{
  return get_dist()->BagImprovement(get_data(), &kFuncEstimate[0], kShrinkage, kDeltaEstimate);
}

//-----------------------------------
// Function: BagData
//
// Returns: none
//
// Description: put data into bags.
//
// Parameters: bool - determines if distribution is pairwise
//    CDistribution ptr - pointer to the distribution + data
//
//-----------------------------------
void CGBMDataContainer::BagData()
{
  get_data().clear_bag();
  get_dist()->BagData(get_data());

}
