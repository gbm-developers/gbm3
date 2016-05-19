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
	static CDistribution* Create(DataDistParams& distParams);

	//---------------------
	// Public destructor
	//---------------------
    virtual ~CPoisson();

    //---------------------
    // Public Functions
    //---------------------
    void ComputeWorkingResponse(const CDataset& data,
    		const double *adF,
				double *adZ);

    double Deviance(const CDataset& data,
    				const double *adF,
                    bool isValidationSet=false);

    double InitF(const CDataset& data);

    void FitBestConstant(const CDataset& data,
    		const double *adF,
			 unsigned long cTermNodes,
			 double* adZ,
			 CTreeComps& treeComps);

    double BagImprovement(const CDataset& data,
			  const double *adF,
			  const double shrinkage,
			  const double* adFadj);

private:
    //----------------------
    // Private Constructors
    //----------------------
    CPoisson();
};

#endif // POISSON_H



