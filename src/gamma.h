//------------------------------------------------------------------------------
//
//  File:       gamma.h
//
//  Description: gamma distribution
//
//------------------------------------------------------------------------------

#ifndef __gamma_h__
#define __gamma_h__

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
	static CDistribution* Create(SEXP radMisc,
								const char* szIRMeasure,
								int& cTrain);

    //---------------------
    // Public destructor
    //---------------------
    virtual ~CGamma();

    //---------------------
    // Public Functions
    //---------------------
    void ComputeWorkingResponse(const CDataset* pData,
    			const double *adF,
				double *adZ);

    double InitF(const CDataset* pData);
    
    void FitBestConstant(const CDataset* pData,
    		const double *adF,
			 unsigned long cTermNodes,
			 double* adZ, CTreeComps* pTreeComps);
    
    double Deviance(const CDataset* pData,
    				const double *adF,
                    bool isValidationSet=false);

    double BagImprovement(const CDataset& data,
    					  const double *adF,
    					  const bag& afInBag,
                          const double shrinkage, const double* adFadj);
private:
    //----------------------
    // Private Constructors
    //----------------------
    CGamma(SEXP radMisc);

};

#endif // __gamma_h__



