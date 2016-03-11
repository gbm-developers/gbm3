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
#include "gbmDataContainer.h"

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
CGBMDataContainer::CGBMDataContainer(SEXP radY, SEXP radOffset, SEXP radX, SEXP raiXOrder,
        SEXP radWeight, SEXP racVarClasses,
        SEXP ralMonotoneVar, SEXP radMisc, const std::string& family, int cTrain, int& cGroups):
        data(radY, radOffset, radX, raiXOrder,
    			radWeight, racVarClasses, ralMonotoneVar, cTrain)
{
	cGroups = -1;

	//Initialize the factory and then use to get the disribution
	DistFactory = new DistributionFactory();

	// Checks for pairwise distribution
	// this should be removed later.
	if(0 == family.compare(0, 8, "pairwise"))
	{
		std::size_t offsetToMeasure = family.find("_");
		if(offsetToMeasure == std::string::npos)
		{
			throw GBM::failure("Unable to locate IR metric required for pairwise");
		}

		const char* szIRMeasure = family.c_str() + offsetToMeasure + 1;
		pDist = DistFactory -> CreateDist("pairwise", radMisc, szIRMeasure, cGroups, cTrain);

	}
	else
	{
		pDist = DistFactory -> CreateDist(family, radMisc, "", cGroups, cTrain);
	}
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
	delete pDist;
	delete DistFactory;
}

//-----------------------------------
// Function: Initialize
//
// Returns: none
//
// Description: Set up the distribution
//
// Parameters: none
//-----------------------------------
void CGBMDataContainer::Initialize()
{
	pDist->Initialize(&data);
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
void CGBMDataContainer::InitializeFunctionEstimate(double& dInitF, unsigned long cLength)
{
	pDist->InitF(&data, dInitF, cLength);
}

//-----------------------------------
// Function: ComputeResiduals
//
// Returns: none
//
// Description: Compute the residuals associated with the distributions loss function.
//
// Parameters: const double ptr - ptr to the function estimates for each predictor
//    CTreeComps ptr - ptr to the tree components container in the gbm.
//
//-----------------------------------
void CGBMDataContainer::ComputeResiduals(const double* adF, CTreeComps* pTreeComp)
{
	pDist->ComputeWorkingResponse(&data, adF,
	                               	pTreeComp->GetGrad(),
	                                pTreeComp->GetBag(),
	                                pTreeComp->GetTrainNo());
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
void CGBMDataContainer::ComputeBestTermNodePreds(const double* adF, CTreeComps* pTreeComp, int& cNodes)
{
	pDist->FitBestConstant(&data, &adF[0],
	                         pTreeComp->GetGrad(),
	                         pTreeComp->GetNodeAssign(),
	                         pTreeComp->GetTrainNo(),
	                         pTreeComp->GetTermNodes(),
	                         (2*cNodes+1)/3, // number of terminal nodes
	                         pTreeComp->GetMinNodeObs(),
	                         pTreeComp->GetBag(),
	                         pTreeComp->GetRespAdj());
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
double CGBMDataContainer::ComputeDeviance(const double* adF, CTreeComps* pTreeComp, bool isValidationSet)
{
	if(!(isValidationSet))
	{
		return pDist->Deviance(&data, adF, pTreeComp->GetTrainNo());
	}
	else
	{
		return pDist->Deviance(&data, adF + pTreeComp->GetTrainNo(), pTreeComp->GetValidNo(), true);
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
double CGBMDataContainer::ComputeBagImprovement(const double* adF, CTreeComps* pTreeComp)
{
	return pDist->BagImprovement(&data, &adF[0],
            pTreeComp->GetRespAdj(),
            pTreeComp->GetBag(),
            pTreeComp->GetLambda(),
            pTreeComp->GetTrainNo());
}

//-----------------------------------
// Function: getDist
//
// Returns: CDistribution ptr
//
// Description: Get pointer to the distribution in use.
//
// Parameters: none
//-----------------------------------
CDistribution* CGBMDataContainer::getDist()
{
	return pDist;
}

//-----------------------------------
// Function: getData
//
// Returns: const Dataset ptr
//
// Description: Get const pointer to dataset.
//
// Parameters: none
//-----------------------------------
const CDataset* CGBMDataContainer::getData()
{
	return &data;
}
