//------------------------------------------------------------------------------
//
//  File:       genericCoxState.h
//
//  Description: abstract class defining the generic methods associated with the CoxPh model.
//
//	Author: 	James Hickey
//------------------------------------------------------------------------------

#ifndef __genericCoxState_h__
#define __genericCoxState_h__

//-----------------------------------
// Definitions
//-----------------------------------
#define frac .00000001
#define recenter 50

//------------------------------
// Includes
//------------------------------
#include "dataset.h"
#include <Rcpp.h>

//------------------------------
// Generic Dispatch Definition
//------------------------------
class GenericCoxState
{
public:
	//----------------------
	// Public Constructors
	//----------------------
	GenericCoxState(){};

	//---------------------
	// Public destructor
	//---------------------
	virtual ~GenericCoxState(){};

	//---------------------
	// Public Functions
	//---------------------
    virtual void ComputeWorkingResponse(const CDataset& data,
    		const double *adF,
				double *adZ)=0;

    virtual void FitBestConstant(const CDataset& data,
    		const double *adF,
			 unsigned long cTermNodes,
			 double* adZ,
			 CTreeComps& treeComps)=0;

    virtual double Deviance(const long cLength, const CDataset& data,
    				const double *adF)=0;

    virtual double BagImprovement(const CDataset& data,
				  const double *adF,
				  const double shrinkage,
				  const double* adFadj)=0;

};
#endif // __genericCoxState_h__
