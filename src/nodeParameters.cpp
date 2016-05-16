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
NodeParams::~NodeParams()
>>>>>>> Merge remote-tracking branch 'pdmOrigin/speed' into trial_stuff
{

}

void NodeParams::ResetSplitProperties(double weightedResiduals, double trainingWeight,
				      unsigned long numObs, double splitValue, unsigned long variableClasses, unsigned long splitVar)

{

  right.weightResid = weightedResiduals;
  right.totalWeight = trainingWeight;
  right.numObs = numObs;

  left.clear();
  missing.clear();

  
  SplitVar = splitVar;
  SplitValue = splitValue;
  ImprovedResiduals = 0.0;
  SplitClass = variableClasses;

  std::fill(adGroupSumZ.begin(), adGroupSumZ.begin() + variableClasses, 0);
  std::fill(adGroupW.begin(), adGroupW.begin() + variableClasses, 0);
  std::fill(acGroupN.begin(), acGroupN.begin() + variableClasses, 0);
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
