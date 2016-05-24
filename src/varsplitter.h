//------------------------------------------------------------------------------
//
//  File:       varsplitter.h
//
//  Description: header for class that splits a node on a particular variable.
//
//------------------------------------------------------------------------------

#ifndef VARSPLITTER_H
#define VARSPLITTER_H

//------------------------------
// Includes
//------------------------------
#include "node.h"
#include "node_parameters.h"
#include <Rcpp.h>

//------------------------------
// Class Definition
//------------------------------
class VarSplitter
{
public:
	//----------------------
	// Public Constructors
	//----------------------
	VarSplitter(unsigned long minNumObs);

	//---------------------
	// Public destructor
	//---------------------
	~VarSplitter();

	//---------------------
	// Public Functions
	//---------------------
	void SetToSplit()
	{
		issplit_ = true;
	};

	 void IncorporateObs(double dX,
			     double dZ,
			     double dW,
			     long lMonotone);

	void Set(CNode& nodeToSplit);
	void ResetForNewVar(unsigned long iWhichVar,
			    long cVarClasses);


	inline double best_improvement() { return bestsplit_.improvement_; }
	inline NodeParams best_split() { return bestsplit_;}
	void SetupNewNodes(CNode& nodeToSplit)
	{
	  nodeToSplit.SplitNode(bestsplit_);
	}

	unsigned long SetAndReturnNumGroupMeans()
	{
	  unsigned long cFiniteMeans = 0;

	  for(unsigned long i=0; i < proposedsplit_.split_class_; i++)
	    {
	      groupMeanAndCat[i].second = i;
	      
	      if(group_weight_[i] != 0.0)
		{
		  groupMeanAndCat[i].first = group_sumresid_[i]/group_weight_[i];
		  cFiniteMeans++;
		}
	      else
		{
		  groupMeanAndCat[i].first = HUGE_VAL;
		}
	    }
	  
	  std::sort(groupMeanAndCat.begin(),
		    groupMeanAndCat.begin() + proposedsplit_.split_class_);

	  return cFiniteMeans;
	}

	void IncrementCategories(unsigned long cat,
				 double predIncrement,
				 double trainWIncrement)
	{
	  group_sumresid_[cat] += predIncrement;
	  group_weight_[cat] += trainWIncrement;
	  group_num_obs_[cat]++;
	}
	
	void UpdateLeftNodeWithCat(long catIndex)
	{

		proposedsplit_.UpdateLeftNode(group_sumresid_[groupMeanAndCat[catIndex].second],
				group_weight_[groupMeanAndCat[catIndex].second],
				group_num_obs_[groupMeanAndCat[catIndex].second]);
	}

	void EvaluateCategoricalSplit();
	void WrapUpCurrentVariable();

	double initial_totalweight;
	double initial_sumresiduals;
	unsigned long initial_numobs;

private:

	//---------------------
	// Private Functions
	//---------------------
	

	//---------------------
	// Private Variables
	//---------------------
	bool issplit_;
	unsigned long min_num_node_obs_;
	double last_xvalue_;
	NodeParams bestsplit_, proposedsplit_;
	std::vector<double> group_sumresid_;
	std::vector<double> group_weight_;
	std::vector<unsigned long> group_num_obs_;

	// Splitting arrays for Categorical variable
	std::vector<std::pair<double, int> > groupMeanAndCat;
};
#endif // VARSPLITTER_H
