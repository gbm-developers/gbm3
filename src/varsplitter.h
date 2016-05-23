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
#include "nodeParameters.h"
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
		fIsSplit = true;
	};

	 void IncorporateObs(double dX,
			     double dZ,
			     double dW,
			     long lMonotone);

	void Set(CNode& nodeToSplit);
	void ResetForNewVar(unsigned long iWhichVar,
			    long cVarClasses);


	inline double best_improvement() { return bestSplit.ImprovedResiduals; }
	inline NodeParams best_split() { return bestSplit;}
	void SetupNewNodes(CNode& nodeToSplit)
	{
	  nodeToSplit.SplitNode(bestSplit);
	}

	unsigned long SetAndReturnNumGroupMeans()
	{
	  unsigned long cFiniteMeans = 0;

	  for(unsigned long i=0; i < proposedSplit.SplitClass; i++)
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
	  
	  std::sort(groupMeanAndCat.begin(),
		    groupMeanAndCat.begin() + proposedSplit.SplitClass);

	  return cFiniteMeans;
	}

	void IncrementCategories(unsigned long cat,
				 double predIncrement,
				 double trainWIncrement)
	{
	  adGroupSumZ[cat] += predIncrement;
	  adGroupW[cat] += trainWIncrement;
	  acGroupN[cat]++;
	}
	
	void UpdateLeftNodeWithCat(long catIndex)
	{

		proposedSplit.UpdateLeftNode(adGroupSumZ[groupMeanAndCat[catIndex].second],
				adGroupW[groupMeanAndCat[catIndex].second],
				acGroupN[groupMeanAndCat[catIndex].second]);
	}

	void EvaluateCategoricalSplit();
	void WrapUpCurrentVariable();

	double dInitTotalW;
	double dInitSumZ;
	unsigned long cInitN;

private:

	//---------------------
	// Private Functions
	//---------------------
	

	//---------------------
	// Private Variables
	//---------------------
	bool fIsSplit;
	unsigned long cMinObsInNode;
	double dLastXValue;
	NodeParams bestSplit, proposedSplit;
	std::vector<double> adGroupSumZ;
	std::vector<double> adGroupW;
	std::vector<unsigned long> acGroupN;

	// Splitting arrays for Categorical variable
	std::vector<std::pair<double, int> > groupMeanAndCat;
};
#endif // VARSPLITTER_H
