//------------------------------------------------------------------------------
//
//  File:       nodeParameters.h
//
//  Description: class implementing the parameters in each node
//
//------------------------------------------------------------------------------
//---------------------
// Includes
//---------------------
#include "nodeParameters.h"
#include "gbmexcept.h"

//---------------------
// Public Functions
//---------------------
NodeParams::NodeParams()
{
	RightNumObs = 0;
	RightTotalWeight = 0.0;
	RightWeightResiduals = 0.0;

	LeftWeightResiduals   = 0.0;
	LeftTotalWeight = 0.0;
	LeftNumObs     = 0;

	MissingWeightResiduals   = 0.0;
	MissingTotalWeight = 0.0;
	MissingNumObs     = 0;

	ImprovedResiduals = 0.0;
	SplitValue = 0.0;
	SplitClass = 0;
	SplitVar = 0;

	adGroupSumZ.resize(1024);
	adGroupW.resize(1024);
	acGroupN.resize(1024);
	groupMeanAndCat.resize(1024);
	aiBestCategory.resize(1024);
}

NodeParams::~NodeParams()
{

}

/*void NodeParams::ResetSplitProperties(double weightedResiduals, double trainingWeight,
									   unsigned long numObs, double splitValue, unsigned long variableClasses, unsigned long splitVar)
{

		RightWeightResiduals   = weightedResiduals;
		RightTotalWeight = trainingWeight;
		RightNumObs     = numObs;


		LeftWeightResiduals   = 0.0;
		LeftTotalWeight = 0.0;
		LeftNumObs     = 0;


		MissingWeightResiduals   = 0.0;
		MissingTotalWeight = 0.0;
		MissingNumObs     = 0;


		SplitVar = splitVar;
		SplitValue = splitValue;
		ImprovedResiduals = 0.0;
		SplitClass = variableClasses;

		//std::cout << variableClasses << "\n";
		std::fill(adGroupSumZ.begin(), adGroupSumZ.begin() + variableClasses, 0);
		std::fill(adGroupW.begin(), adGroupW.begin() + variableClasses, 0);
		std::fill(acGroupN.begin(), acGroupN.begin() + variableClasses, 0);

}


void NodeParams::UpdateMissingNode(double predIncrement, double trainWIncrement, long numIncrement)
{
	// Move data point from right node to missing
	MissingWeightResiduals += predIncrement;
	MissingTotalWeight += trainWIncrement;
	MissingNumObs += numIncrement;

	RightWeightResiduals -= predIncrement;
	RightTotalWeight -= trainWIncrement;
	RightNumObs -= numIncrement;
}

void NodeParams::UpdateLeftNode(double predIncrement, double trainWIncrement, long numIncrement)
{
	// Move data point from right node to left node
	LeftWeightResiduals += predIncrement;
	LeftTotalWeight += trainWIncrement;
	LeftNumObs += numIncrement;

	RightWeightResiduals -= predIncrement;
	RightTotalWeight -= trainWIncrement;
	RightNumObs -= numIncrement;

}

void NodeParams::NodeGradResiduals()
{
	// Returns weighted

	double dTemp = 0.0;
	double dResult = 0.0;

	// Only need to look at left and right
	if(MissingNumObs == 0.0)
	{
		dTemp = LeftWeightResiduals/LeftTotalWeight - RightWeightResiduals/RightTotalWeight;
		dResult = LeftTotalWeight*RightTotalWeight*dTemp*dTemp/(LeftTotalWeight+RightTotalWeight);
	}
	else
	{
		// Grad - left/right
		dTemp = LeftWeightResiduals/LeftTotalWeight - RightWeightResiduals/RightTotalWeight;
		dResult += LeftTotalWeight*RightTotalWeight*dTemp*dTemp;

		// Grad - left/missing
		dTemp = LeftWeightResiduals/LeftTotalWeight - MissingWeightResiduals/MissingTotalWeight;
		dResult += LeftTotalWeight*MissingTotalWeight*dTemp*dTemp;

		// Grad - right/missing
		dTemp = RightWeightResiduals/RightTotalWeight - MissingWeightResiduals/MissingTotalWeight;
		dResult += RightTotalWeight*MissingTotalWeight*dTemp*dTemp;
		dResult /= (LeftTotalWeight + RightTotalWeight + MissingTotalWeight);
	}

	// Update current residuals
	ImprovedResiduals = dResult;
}

bool NodeParams::SplitIsCorrMonotonic(long specifyMonotone)
{
	double weightedGrad = RightWeightResiduals * LeftTotalWeight- LeftWeightResiduals * RightTotalWeight;
	return (specifyMonotone == 0 || specifyMonotone * weightedGrad > 0);
}

bool NodeParams::HasMinNumOfObs(long minObsInNode)
{
	return ((LeftNumObs >= minObsInNode) &&
				(RightNumObs >= minObsInNode));
}

void NodeParams::IncrementCategories(unsigned long cat, double predIncrement, double trainWIncrement)
{
	adGroupSumZ[cat] += predIncrement;
	adGroupW[cat] += trainWIncrement;
	acGroupN[cat]++;
}


unsigned long NodeParams::SetAndReturnNumGroupMeans()
{
	unsigned long cFiniteMeans = 0;

	for(long i=0; i < SplitClass; i++)
	{
	  groupMeanAndCat[i].second = i;

	  if(adGroupW[i] != 0.0)
	  {
		  groupMeanAndCat[i].first = adGroupSumZ[i]/adGroupW[i];
		  cFiniteMeans++;
	  }
	  else
	  {
		  groupMeanAndCat[i].first = HUGE_VAL;
	  }
	}

  std::sort(groupMeanAndCat.begin(), groupMeanAndCat.begin() + SplitClass);

  return cFiniteMeans;
}

void NodeParams::UpdateLeftNodeWithCat(long catIndex)
{

	UpdateLeftNode(adGroupSumZ[groupMeanAndCat[catIndex].second],
			adGroupW[groupMeanAndCat[catIndex].second],
			acGroupN[groupMeanAndCat[catIndex].second]);
}*/
