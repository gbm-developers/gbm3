//------------------------------------------------------------------------------
//
//  File:       tweedie.h
//
//  Description: tweedie distribution for .
//
//------------------------------------------------------------------------------

#ifndef __tweedie_h__
#define __tweedie_h__

//------------------------------
// Includes
//------------------------------
#include "distribution.h"
#include <Rmath.h>
#include <Rcpp.h>
#include <memory>

//------------------------------
// Class definition
//------------------------------
class CTweedie : public CDistribution
{

public:
	//---------------------
	// Factory Function
	//---------------------
	static CDistribution* Create(SEXP radMisc,
								const char* szIRMeasure, int& cTrain);

	//---------------------
	// Public destructor
	//---------------------
    virtual ~CTweedie();

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
			 double* adZ,
			 CTreeComps* pTreeComps);

    double Deviance(const CDataset* pData,
    				const double *adF,
                    bool isValidationSet=false);
    
    double BagImprovement(const CDataset& data,
    					  const double *adF,
    					  const bag& afInBag,
                          const double shrinkage,
                          const double* adFadj);
private:
    //----------------------
    // Private Constructors
    //----------------------
    CTweedie(SEXP radMisc);

	//-------------------
	// Private Variables
	//-------------------
    double dPower;
};

#endif // _tweedie_h__

