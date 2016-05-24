//------------------------------------------------------------------------------
//
//  File:       genericNodeStrategy.h
//
//  Description: abstract class defining the generic node strategy methods.
//
//	Author: 	James Hickey
//------------------------------------------------------------------------------

#ifndef GENERICNODESTRATEGY_H
#define GENERICNODESTRATEGY_H

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
	virtual void Adjust(unsigned long min_num_node_obs)=0;
	virtual void Predict(const CDataset &kData,
		    unsigned long rownum,
		    double &delta_estimate)=0;
	virtual void GetVarRelativeInfluence(double* relative_influence)=0;
	virtual void PrintSubTree(unsigned long indent)=0;
	virtual signed char WhichNode(const CDataset& kData, unsigned long obs_num)=0;
	virtual void TransferTreeToRList
	(
	    int &nodeid,
	    const CDataset &kData,
	    int* splitvar,
	    double* splitvalues,
	    int* leftnodes,
	    int* rightnodes,
	    int* missingnodes,
	    double* error_reduction,
	    double* weights,
	    double* predictions,
	    VEC_VEC_CATEGORIES &splitcodes_vec,
	    int prev_categorical_splits,
	    double shrinkage
	)=0;

};
#endif // GENERICNODESTRATEGY_H
