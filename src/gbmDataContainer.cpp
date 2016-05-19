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
CGBMDataContainer::CGBMDataContainer(DataDistParams dataDistConfig):
        data(dataDistConfig)
{

	//Initialize the factory and then use to get the disribution
	DistFactory = new DistributionFactory();
	pDist = DistFactory -> CreateDist(dataDistConfig);

	// Set up multi map
	for(long i = 0; i < dataDistConfig.patId.size(); i++)
	{
		patIdToRow.insert(pair<int, int>(dataDistConfig.patId(i), i));
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
// Returns: nonei = 0; i < data.getNumUniquePatient()
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
		pair<multimap<int, int>::iterator, multimap<int, int>::iterator> keyRange;
		multimap<int, int>::iterator patIt, rowIt;

		long finalRowBeforeStopBagging = 0;
		bool bagVal;

		// Bag via patient id  - loop over patients
		for(patIt = patIdToRow.begin();
				patIt != patIdToRow.end();
				patIt = rowIt)
		{

			// Check if we've filled the bag or have left the training set
			if((i >= data.GetNumPatientsInTraining()) || (cBagged >= data.GetTotalInBag())) break;

			keyRange = patIdToRow.equal_range(patIt->first);

			// Check if that patient should be bagged
			if(unif_rand() * (data.GetNumPatientsInTraining()-i) < data.GetTotalInBag() - cBagged)
			{
				bagVal = true;
				cBagged++;
			}
			else
			{
				bagVal = false;
			}

			// Loop over rows and set to bool value
			for(rowIt = keyRange.first; rowIt != keyRange.second; ++rowIt)
			{
				finalRowBeforeStopBagging = (*rowIt).second;
				data.SetBagElem((*rowIt).second, bagVal);
			}

			// Increment patient number
			i += 1;
		}

		data.FillRemainderOfBag(finalRowBeforeStopBagging);

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

			const double dGroup = static_cast<CPairwise*>(pDist)->adGroup[i];
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
