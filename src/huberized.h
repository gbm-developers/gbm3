//------------------------------------------------------------------------------
//
//  File:       huberized.h
//
//  Description:   huberized hinge loss object for GBM.
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//------------------------------------------------------------------------------

#ifndef __huberized_h__
#define __huberized_h__

//------------------------------
// Includes
//------------------------------
#include "distribution.h"
#include "buildinfo.h"
#include <memory>

//------------------------------
// Class definition
//------------------------------
class CHuberized : public CDistribution
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
    virtual ~CHuberized();

    //---------------------
    // Public Functions
    //---------------------
    void ComputeWorkingResponse(const CDataset* pData,
    			const double *adF,
				double *adZ);

    double Deviance(const CDataset* pData,
    				const double *adF,
                    bool isValidationSet=false);

    double InitF(const CDataset* pData);

    void FitBestConstant(const CDataset* pData,
    		const double *adF,
			 unsigned long cTermNodes,
			 double* adZ,
			 CTreeComps* pTreeComps);

    double BagImprovement(const CDataset& data,
    					  const double *adF,
    					  const bag& afInBag,
                          const double shrinkage,
                          const double* adFadj);

private:
    //----------------------
    // Private Constructors
    //----------------------
    CHuberized(SEXP radMisc);
};

#endif // __huberized_h__



