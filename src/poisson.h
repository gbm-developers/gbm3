//------------------------------------------------------------------------------
//
//  File:       poisson.h
//
//  Description:   poisson object for GBM.
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef POISSON_H
#define POISSON_H

//------------------------------
// Includes
//------------------------------
#include "distribution.h"
#include <Rmath.h>
#include <memory>

//------------------------------
// Class definition
//------------------------------
class CPoisson : public CDistribution
{

public:
	//---------------------
	// Factory Function
	//---------------------
	static CDistribution* Create(DataDistParams& distparams);

	//---------------------
	// Public destructor
	//---------------------
    virtual ~CPoisson();

    //---------------------
    // Public Functions
    //---------------------
    void ComputeWorkingResponse(const CDataset& kData,
    		const double *kFuncEstimate,
				double *residual);

    double Deviance(const CDataset&kData,
    				const double *kFuncEstimate,
                    bool is_validationset=false);

    double InitF(const CDataset&kData);

    void FitBestConstant(const CDataset&kData,
    		const double *kFuncEstimate,
			 unsigned long num_terminalnodes,
			 double* residual,
			 CTreeComps& treecomps);

    double BagImprovement(const CDataset&kData,
			  const double *kFuncEstimate,
			  const double kShrinkage,
			  const double* kFuncEstimateadj);

private:
    //----------------------
    // Private Constructors
    //----------------------
    CPoisson();
};

#endif // POISSON_H



