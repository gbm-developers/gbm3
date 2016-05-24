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
  NodeDef() : numobs_(0), weightresid_(0), totalweight_(0) {};
  
  NodeDef(double weightResid, double totalWeight, unsigned long numObs) :
  numobs_(numObs), weightresid_(weightResid), totalweight_(totalWeight) {};
  
  //---------------------
  // Public Functions
  //---------------------
  void clear()
  {
    numobs_ = 0;
    weightresid_ = totalweight_ = 0;
  };
  
  void increment(const double pred, const double trainWeight, long num)
  {
    weightresid_ += pred;
    totalweight_ += trainWeight;
    numobs_ += num;
  };
  
  double prediction() const
  {
    return weightresid_ / totalweight_;
  };
  
  double unweighted_gradient(const NodeDef& other) const
  {
	return weightresid_ * other.totalweight_ - other.weightresid_ * totalweight_;
  };

  double variance_reduction(const NodeDef& other) const
  {
	const double predictionDiff = prediction() - other.prediction();
    return totalweight_* other.totalweight_*predictionDiff*predictionDiff;
  };

  bool has_min_obs(unsigned long minObsInNode) const
  {
	return (numobs_ >= minObsInNode);
  }

  bool has_obs() const
  {
    return numobs_;
  }

  double sum_weights(const NodeDef& a) const
  {
    return totalweight_ + a.totalweight_;
  }
  
  double sum_weights(const NodeDef& a, const NodeDef& b) const
  {
    return totalweight_ + a.sum_weights(b);
  }

  double get_totalweight() const
  {
    return totalweight_;
  }

  long get_num_obs() const
  {
    return numobs_;
  }
  
private:
  
  unsigned long numobs_;
  double weightresid_;
  double totalweight_;

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
	NodeParams() : category_ordering_(1024) {};

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
	  missing_.increment(predIncrement, trainWIncrement, numIncrement);
	  right_.increment(-predIncrement, -trainWIncrement, -numIncrement);
	}
	void UpdateLeftNode(double predIncrement, double trainWIncrement, long numIncrement = 1)
	{
		// Move data point from right node to left node
		left_.increment(predIncrement, trainWIncrement, numIncrement);
		right_.increment(-predIncrement, -trainWIncrement, -numIncrement);
	}
	inline double get_improvement() { return improvement_;};
	bool split_is_correct_monotonicity(long specifyMonotone)
	{
		return (
			(specifyMonotone == 0) ||
			((specifyMonotone * right_.unweighted_gradient(left_)) > 0)
			);
	}
	void NodeGradResiduals()
	{
	  // Only need to look at left and right
	  if(!missing_.has_obs())
	    {
	      improvement_ = left_.variance_reduction(right_) /
		left_.sum_weights(right_);
	    }
	  else
	    {
	      // Grad - left/right
	      improvement_ =
		(left_.variance_reduction(right_) +
		 left_.variance_reduction(missing_) +
		 right_.variance_reduction(missing_)) /
		left_.sum_weights(right_, missing_);
	    }
	};
	
	bool has_min_num_obs(unsigned long minObsInNode)
	{
		return (left_.has_min_obs(minObsInNode) &&
			  right_.has_min_obs(minObsInNode));
	}
	inline void SetBestCategory(std::vector<std::pair<double, int> >& groupMeanAndCat)
	{

	  int count = 0;
	  category_ordering_.resize(groupMeanAndCat.size());
	  for(std::vector<std::pair<double, int> >::const_iterator it = groupMeanAndCat.begin();
	      it != groupMeanAndCat.end();
	      ++it)
	    {
	      category_ordering_[count] = it->second;
	      count++;
	    }
	};

	bool has_missing() const
	{
	  return missing_.has_obs();
	};
	//---------------------
	// Public Variables
	//---------------------
	// Left Node Definition
	NodeDef left_, right_, missing_;

	// Splitting values
	double split_value_; // Continuous Split Value
	unsigned long split_var_; // Which feature to split on
	unsigned long split_class_; // Categorical Split Value
	std::vector<int> category_ordering_; // Vector of levels ordering
	double improvement_;
};

#endif // NODEPARAMETERS_H
