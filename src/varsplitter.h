//------------------------------------------------------------------------------
//
//  File:       varsplitter.h
//
//  Description: header for class that splits a node on a particular variable.
//
//------------------------------------------------------------------------------

#ifndef __varsplitter_h__
#define __varsplitter_h__

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


	inline double BestImprovement() { return bestSplit.ImprovedResiduals; }
	inline NodeParams GetBestSplit() { return bestSplit;}
	void SetupNewNodes(CNode& nodeToSplit)
	{
		nodeToSplit.SplitNode(bestSplit);
	}

	void EvaluateCategoricalSplit();
	void WrapUpCurrentVariable();

	double dInitTotalW;
	double dInitSumZ;
	unsigned long cInitN;


private:

	unsigned long cMinObsInNode;

	double dLastXValue;

	NodeParams bestSplit, proposedSplit;

	//---------------------
	// Private Functions
	//---------------------
	

	//---------------------
	// Private Variables
	//---------------------
	bool fIsSplit;



};
#endif // __varplitter_h__
