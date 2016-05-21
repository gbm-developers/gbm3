//------------------------------------------------------------------------------
//
//  File:       gamma.h
//
//  Description: gamma distribution
//
//------------------------------------------------------------------------------

#ifndef GAMMA_H
#define GAMMA_H

//------------------------------
// Includes
//------------------------------
#include "distribution.h"
#include <Rmath.h>
#include <memory>

//------------------------------
// Class definition
//------------------------------
class CGamma : public CDistribution
{

public:

	//---------------------
	// Factory Function
	//---------------------
	static CDistribution* Create(DataDistParams& distParams);

    //---------------------
    // Public destructor
    //---------------------
    virtual ~CGamma();

    //---------------------
    // Public Functions
    //---------------------
    void ComputeWorkingResponse(const CDataset& data,
    			const double *adF,
				double *adZ);

    double InitF(const CDataset& data);
    
    void FitBestConstant(const CDataset& data,
    		const double *adF,
			 unsigned long cTermNodes,
			 double* adZ, CTreeComps& treeComps);
    
    double Deviance(const CDataset& data,
    				const double *adF,
                    bool isValidationSet=false);

    double BagImprovement(const CDataset& data,
			  const double *adF,
                          const double shrinkage, const double* adFadj);
private:
    //----------------------
    // Private Constructors
    //----------------------
    CGamma();

};

#endif // GAMMA_H



