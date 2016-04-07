//------------------------------------------------------------------------------
//
//  File:       genericNodeStrategy.h
//
//  Description: abstract class defining the generic node strategy methods.
//
//	Author: 	James Hickey
//------------------------------------------------------------------------------

#ifndef __genericNodeStrategy_h__
#define __genericNodeStrategy_h__

//------------------------------
// Includes
//------------------------------
#include "dataset.h"
#include <Rcpp.h>



//------------------------------
// Generic Dispatch Definition
//------------------------------
class GenericNodeStrategy
{
public:
	//----------------------
	// Public Constructors
	//----------------------
	GenericNodeStrategy(){};

	//---------------------
	// Public destructor
	//---------------------
	virtual ~GenericNodeStrategy(){};

	//---------------------
	// Public Functions
	//---------------------
	virtual void Adjust(unsigned long cMinObsInNode)=0;
	virtual void Predict(const CDataset &data,
		    unsigned long iRow,
		    double &dFadj)=0;
	virtual void GetVarRelativeInfluence(double* adRelInf)=0;
	virtual void PrintSubTree(unsigned long Indent)=0;
	virtual signed char WhichNode(const CDataset& data, unsigned long iObs)=0;
	virtual void TransferTreeToRList
	(
	    int &iNodeID,
	    const CDataset &data,
	    int *aiSplitVar,
	    double *adSplitPoint,
	    int *aiLeftNode,
	    int *aiRightNode,
	    int *aiMissingNode,
	    double *adErrorReduction,
	    double *adWeight,
	    double *adPred,
	    VEC_VEC_CATEGORIES &vecSplitCodes,
	    int cCatSplitsOld,
	    double dShrinkage
	)=0;

};
#endif // __genericNodeStrategy_h__
