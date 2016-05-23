//------------------------------------------------------------------------------
//
//  File:       nodeParameters.h
//
//  Description:  header  contains the parameters used to split a node
//
//------------------------------------------------------------------------------

#ifndef NODEPARAMETERS_H
#define NODEPARAMETERS_H

//------------------------------
// Includes
//------------------------------
#include <Rcpp.h>

//------------------------------
// Struct Definition
//------------------------------
struct NodeDef
{

  //----------------------
  // Public Constructors
  //----------------------
  NodeDef() : numObs(0), weightResid(0), totalWeight(0) {};
  
  NodeDef(double weightResid, double totalWeight, unsigned long numObs) :
  numObs(numObs), weightResid(weightResid), totalWeight(totalWeight) {};
  
  //---------------------
  // Public Functions
  //---------------------
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
  
  double unweighted_gradient(const NodeDef& other) const
  {
	return weightResid * other.totalWeight - other.weightResid * totalWeight;
  };

  double variance_reduction(const NodeDef& other) const
  {
	const double predictionDiff = prediction() - other.prediction();
    return totalWeight* other.totalWeight*predictionDiff*predictionDiff;
  };

  bool has_min_obs(unsigned long minObsInNode) const
  {
	return (numObs >= minObsInNode);
  }

  bool has_obs() const
  {
    return numObs;
  }

  double sum_weights(const NodeDef& a) const
  {
    return totalWeight + a.totalWeight;
  }
  
  double sum_weights(const NodeDef& a, const NodeDef& b) const
  {
    return totalWeight + a.sum_weights(b);
  }

  double get_totalweight() const
  {
    return totalWeight;
  }

  long get_num_obs() const
  {
    return numObs;
  }
  
private:
  
  unsigned long numObs;
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
	NodeParams() : aiBestCategory(1024) {};

	//---------------------
	// Public destructor
	//---------------------
        ~NodeParams();

	//---------------------
	// Public Functions
	//---------------------
	void ResetSplitProperties(double weightedResiduals, double trainingWeight, unsigned long numObs,
				  double splitValue = -HUGE_VAL, unsigned long variableClasses=1, unsigned long splitVar = UINT_MAX);
	void UpdateMissingNode(double predIncrement, double trainWIncrement, long numIncrement = 1)
	{
	  // Move data point from right node to missing
	  missing.increment(predIncrement, trainWIncrement, numIncrement);
	  right.increment(-predIncrement, -trainWIncrement, -numIncrement);
	}
	void UpdateLeftNode(double predIncrement, double trainWIncrement, long numIncrement = 1)
	{
		// Move data point from right node to left node
		left.increment(predIncrement, trainWIncrement, numIncrement);
		right.increment(-predIncrement, -trainWIncrement, -numIncrement);
	}
	inline double get_improvement() { return ImprovedResiduals;};
	bool split_is_correct_monotonicity(long specifyMonotone)
	{
		return (
			(specifyMonotone == 0) ||
			((specifyMonotone * right.unweighted_gradient(left)) > 0)
			);
	}
	void NodeGradResiduals()
	{
	  // Only need to look at left and right
	  if(!missing.has_obs())
	    {
	      ImprovedResiduals = left.variance_reduction(right) /
		left.sum_weights(right);
	    }
	  else
	    {
	      // Grad - left/right
	      ImprovedResiduals =
		(left.variance_reduction(right) +
		 left.variance_reduction(missing) +
		 right.variance_reduction(missing)) /
		left.sum_weights(right, missing);
	    }
	};
	
	bool has_min_num_obs(unsigned long minObsInNode)
	{
		return (left.has_min_obs(minObsInNode) &&
			  right.has_min_obs(minObsInNode));
	}
	inline void setBestCategory(std::vector<std::pair<double, int> >& groupMeanAndCat)
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

	bool has_missing() const
	{
	  return missing.has_obs();
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
};

#endif // NODEPARAMETERS_H
