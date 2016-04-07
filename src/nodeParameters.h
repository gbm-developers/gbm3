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
	double GetImprovement() { return ImprovedResiduals;};
	bool SplitIsCorrMonotonic(long specifyMonotone);
	void NodeGradResiduals();
	bool HasMinNumOfObs(long minObsInNode);
	void setBestCategory()
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
	NodeParams& operator=(const NodeParams rhs)
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

#endif // __nodeParameters_h__
