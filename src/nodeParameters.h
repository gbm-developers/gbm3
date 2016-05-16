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

struct NodeDef
{
  NodeDef() : numObs(0), weightResid(0), totalWeight(0) {};

  NodeDef(double weightResid, double totalWeight, long numObs) :
    weightResid(weightResid), totalWeight(totalWeight), numObs(numObs) {};
  
  void clear()
  {
    numObs = 0;
    weightResid = totalWeight = 0;
  };

  void increment(const double pred, const double trainWeight, long num)
  {
    weightResid += pred;
    totalWeight += trainWeight;
    numObs += num;
  };

  double prediction() const
  {
    return weightResid / totalWeight;
  };
  
  double unweightedGradient(const NodeDef& other) const
  {
    const double tmp = prediction() - other.prediction();
    return totalWeight * other.totalWeight * tmp * tmp;
  };

  bool hasMinObs(long minObsInNode) const
  {
    return (numObs >= minObsInNode);
  }
  
  long numObs;
  double weightResid;
  double totalWeight;

};

//------------------------------
// Class Definition
//------------------------------
class NodeParams
{
public:
	//----------------------
	// Public Constructors
	//----------------------
	NodeParams() : adGroupSumZ(1024), adGroupW(1024), acGroupN(1024), groupMeanAndCat(1024), aiBestCategory(1024) {};

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

	NodeParams& operator=(const NodeParams& rhs)
	{
	  right = rhs.right;
	  left = rhs.left;
	  missing = rhs.missing;

	  SplitValue = rhs.SplitValue;
	  SplitVar = rhs.SplitVar;
	  SplitClass = rhs.SplitClass;
	  ImprovedResiduals = rhs.ImprovedResiduals;
	  
	  // Copy best category
	  aiBestCategory = rhs.aiBestCategory;
	  return *this;
	}
	bool hasMissing() const
	{
		return missing.numObs >= 0;
	};
	//---------------------
	// Public Variables
	//---------------------
	// Left Node Definition
	NodeDef left, right, missing;

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
