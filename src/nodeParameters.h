//------------------------------------------------------------------------------
//
//  File:       nodeParameters.h
//
//  Description:  header  contains the parameters used to split a node
//
//------------------------------------------------------------------------------

#ifndef __nodeParameters_h__
#define __nodeParameters_h__

//------------------------------
// Includes
//------------------------------
#include <Rcpp.h>

//------------------------------
// Class Definition
//------------------------------
class NodeParams
{
public:
	//----------------------
	// Public Constructors
	//----------------------
	NodeParams();

	//---------------------
	// Public destructor
	//---------------------
    ~NodeParams();

	//---------------------
	// Public Functions
	//---------------------
	void ResetSplitProperties(double weightedResiduals, double trainingWeight, unsigned long numObs,
							 double splitValue = -HUGE_VAL, unsigned long variableClasses=1, unsigned long splitVar = UINT_MAX);
	void UpdateMissingNode(double predIncrement, double trainWIncrement, long numIncrement = 1);
	void UpdateLeftNode(double predIncrement, double trainWIncrement, long numIncrement = 1);
	void UpdateLeftNodeWithCat(long catIndex);
	void IncrementCategories(unsigned long cat, double predIncrement, double trainWIncrement);
	unsigned long SetAndReturnNumGroupMeans();
	inline double GetImprovement() { return ImprovedResiduals;};
	bool SplitIsCorrMonotonic(long specifyMonotone);
	void NodeGradResiduals();
	bool HasMinNumOfObs(long minObsInNode);
	inline void setBestCategory()
	{
		int count = 0;
		aiBestCategory.resize(groupMeanAndCat.size());
		for(std::vector<std::pair<double, int> >::const_iterator it = groupMeanAndCat.begin();
				it != groupMeanAndCat.end();
				++it)
		{
			aiBestCategory[count] = it->second;
			count++;
		}
	};
	inline NodeParams& operator=(const NodeParams rhs)
	{
		RightWeightResiduals = rhs.RightWeightResiduals;
		RightTotalWeight = rhs.RightTotalWeight;
		RightNumObs = rhs.RightNumObs;

		LeftWeightResiduals = rhs.LeftWeightResiduals;
		LeftTotalWeight = rhs.LeftTotalWeight;
		LeftNumObs = rhs.LeftNumObs;

		MissingWeightResiduals = rhs.MissingWeightResiduals;
		MissingTotalWeight = rhs.MissingTotalWeight;
		MissingNumObs = rhs.MissingNumObs;

		SplitValue = rhs.SplitValue;
		SplitVar = rhs.SplitVar;
		SplitClass = rhs.SplitClass;
		ImprovedResiduals = rhs.ImprovedResiduals;

		// Copy best category
		aiBestCategory.resize(rhs.aiBestCategory.size(), 0);
		std::copy(rhs.aiBestCategory.begin(), rhs.aiBestCategory.end(), aiBestCategory.begin());
		return *this;
	}
	//---------------------
	// Public Variables
	//---------------------
	// Left Node Definition
	double LeftWeightResiduals;
	double LeftTotalWeight;
	long LeftNumObs;

	// Right Node Definition
	double RightWeightResiduals;
	double RightTotalWeight;
	long RightNumObs;

	// Missing Node Definition
	double MissingWeightResiduals;
	double MissingTotalWeight;
	long MissingNumObs;

	// Splitting values
	double SplitValue; // Continuous Split Value
	unsigned long SplitVar; // Which feature to split on
	unsigned long SplitClass; // Categorical Split Value
    std::vector<int> aiBestCategory; // Vector of levels ordering
	double ImprovedResiduals;

	// Splitting arrays for Categorical variable
	std::vector<double> adGroupSumZ;
	std::vector<double> adGroupW;
	std::vector<unsigned long> acGroupN;
	std::vector<std::pair<double, int> > groupMeanAndCat;

};

inline void NodeParams::ResetSplitProperties(double weightedResiduals, double trainingWeight,
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


inline void NodeParams::UpdateMissingNode(double predIncrement, double trainWIncrement, long numIncrement)
{
	// Move data point from right node to missing
	MissingWeightResiduals += predIncrement;
	MissingTotalWeight += trainWIncrement;
	MissingNumObs += numIncrement;

	RightWeightResiduals -= predIncrement;
	RightTotalWeight -= trainWIncrement;
	RightNumObs -= numIncrement;
}

inline void NodeParams::UpdateLeftNode(double predIncrement, double trainWIncrement, long numIncrement)
{
	// Move data point from right node to left node
	LeftWeightResiduals += predIncrement;
	LeftTotalWeight += trainWIncrement;
	LeftNumObs += numIncrement;

	RightWeightResiduals -= predIncrement;
	RightTotalWeight -= trainWIncrement;
	RightNumObs -= numIncrement;

}

inline void NodeParams::NodeGradResiduals()
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

inline bool NodeParams::SplitIsCorrMonotonic(long specifyMonotone)
{
	double weightedGrad = RightWeightResiduals * LeftTotalWeight- LeftWeightResiduals * RightTotalWeight;
	return (specifyMonotone == 0 || specifyMonotone * weightedGrad > 0);
}

inline bool NodeParams::HasMinNumOfObs(long minObsInNode)
{
	return ((LeftNumObs >= minObsInNode) &&
				(RightNumObs >= minObsInNode));
}

inline void NodeParams::IncrementCategories(unsigned long cat, double predIncrement, double trainWIncrement)
{
	adGroupSumZ[cat] += predIncrement;
	adGroupW[cat] += trainWIncrement;
	acGroupN[cat]++;
}


inline unsigned long NodeParams::SetAndReturnNumGroupMeans()
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

inline void NodeParams::UpdateLeftNodeWithCat(long catIndex)
{

	UpdateLeftNode(adGroupSumZ[groupMeanAndCat[catIndex].second],
			adGroupW[groupMeanAndCat[catIndex].second],
			acGroupN[groupMeanAndCat[catIndex].second]);
}

#endif // __nodeParameters_h__
