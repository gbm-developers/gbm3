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
CGBMDataContainer::CGBMDataContainer(DataDistParams dataDistConfig):
        data(dataDistConfig)
{

	//Initialize the factory and then use to get the disribution
	DistFactory = new DistributionFactory();

	// Checks for pairwise distribution
	// this should be removed later.
	if(0 == dataDistConfig.family.compare(0, 8, "pairwise"))
	{
		std::size_t offsetToMeasure = dataDistConfig.family.find("_");
		if(offsetToMeasure == std::string::npos)
		{
			throw GBM::failure("Unable to locate IR metric required for pairwise");
		}

		const char* szIRMeasure = dataDistConfig.family.c_str() + offsetToMeasure + 1;
		pDist = DistFactory -> CreateDist("pairwise", dataDistConfig.misc, szIRMeasure, dataDistConfig.cTrain);

	}
	else
	{
		pDist = DistFactory -> CreateDist(dataDistConfig.family, dataDistConfig.misc, "", dataDistConfig.cTrain);
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
double CGBMDataContainer::InitialFunctionEstimate()
{
	return pDist->InitF(&data);
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
void CGBMDataContainer::ComputeResiduals(const double* adF, double* adZ)
{
	pDist->ComputeWorkingResponse(&data, adF, adZ);
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
void CGBMDataContainer::ComputeBestTermNodePreds(const double* adF, double* adZ, CTreeComps* pTreeComp)
{
	pDist->FitBestConstant(&data, &adF[0],
	                         (2*pTreeComp->GetSizeOfTree()+1)/3, // number of terminal nodes
	                         &adZ[0], pTreeComp);
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
double CGBMDataContainer::ComputeDeviance(const double* adF, bool isValidationSet)
{
	if(!(isValidationSet))
	{
		return pDist->Deviance(&data, adF);
	}
	else
	{
		return pDist->Deviance(&data, adF + data.get_trainSize(), true);
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
double CGBMDataContainer::ComputeBagImprovement(const double* adF, const double shrinkage, const double* adFadj)
{
	return pDist->BagImprovement(data, &adF[0], data.GetBag(),  shrinkage, adFadj);
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
CDataset* CGBMDataContainer::getData()
{
	return &data;
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
	unsigned long i = 0;
	unsigned long cBagged = 0;

	// randomly assign observations to the Bag
	if (pDist->GetNumGroups() < 0)
	{
		// regular instance based training
		for(i=0; i<data.get_trainSize() && (cBagged < data.GetTotalInBag()); i++)
		{
			if(unif_rand() * (data.get_trainSize()-i) < data.GetTotalInBag() - cBagged)
			{
				data.SetBagElem(i, true);
				cBagged++;
			}
			else
			{
				 data.SetBagElem(i, false);
			}
		}

		data.FillRemainderOfBag(i);
	}
	else
	{
		// for pairwise training, sampling is per group
		// therefore, we will not have exactly cTotalInBag instances

		double dLastGroup = -1;
		bool fChosen = false;
		unsigned int cBaggedGroups = 0;
		unsigned int cSeenGroups   = 0;
		unsigned int cTotalGroupsInBag = (unsigned long)(data.GetBagFraction() * pDist->GetNumGroups());
		if (cTotalGroupsInBag <= 0)
		{
			cTotalGroupsInBag = 1;
		}

		for(i=0; i< data.get_trainSize(); i++)
		{

			const double dGroup = pDist->misc_ptr(true)[i];

			if(dGroup != dLastGroup)
			{
				if (cBaggedGroups >= cTotalGroupsInBag)
				{
					break;
				}

				// Group changed, make a new decision
				fChosen = (unif_rand()*(pDist->GetNumGroups() - cSeenGroups) <
			   cTotalGroupsInBag - cBaggedGroups);
				if(fChosen)
				{
					cBaggedGroups++;
				}
				dLastGroup = dGroup;
				cSeenGroups++;
			}
			if(fChosen)
			{
				data.SetBagElem(i, true);
				cBagged++;
			}
			else
			{
				data.SetBagElem(i, false);
			}
		}

		// the remainder is not in the bag
		data.FillRemainderOfBag(i);
	}
}
