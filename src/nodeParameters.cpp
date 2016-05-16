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


void NodeParams::UpdateMissingNode(double predIncrement, double trainWIncrement, long numIncrement)
{
  // Move data point from right node to missing
  missing.increment(predIncrement, trainWIncrement, numIncrement);
  right.increment(-predIncrement, -trainWIncrement, -numIncrement);
}

void NodeParams::UpdateLeftNode(double predIncrement, double trainWIncrement, long numIncrement)
{
  // Move data point from right node to left node
  left.increment(predIncrement, trainWIncrement, numIncrement);
  right.increment(-predIncrement, -trainWIncrement, -numIncrement);
}

void NodeParams::NodeGradResiduals()
{
  // Returns weighted
  
  // Only need to look at left and right
  if(missing.numObs == 0)
    {
      ImprovedResiduals = left.unweightedGradient(right) /
	(left.totalWeight + right.totalWeight);
    }
  else
    {
      // Grad - left/right
      
      ImprovedResiduals =
	(left.unweightedGradient(right) +
	 left.unweightedGradient(missing) +
	 right.unweightedGradient(missing)) /
	(left.totalWeight + right.totalWeight + missing.totalWeight);
    }
}

bool NodeParams::SplitIsCorrMonotonic(long specifyMonotone)
{
  double weightedGrad = right.weightResid * left.totalWeight -
    left.weightResid * right.totalWeight;
  return (specifyMonotone == 0 || specifyMonotone * weightedGrad > 0);
}

bool NodeParams::HasMinNumOfObs(long minObsInNode)
{
  return (left.hasMinObs(minObsInNode) &&
	  right.hasMinObs(minObsInNode));
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
