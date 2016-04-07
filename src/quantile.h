//------------------------------------------------------------------------------
//  File:       quantile.h
//
//  Contents:   quantile regression for GBM.
//
//  History:    10/8/2006   Created by Brian Kriegler (bk@stat.ucla.edu)
//              6/11/2007   gregr merged with official gbm
//
//------------------------------------------------------------------------------

#ifndef __quantile_h__
#define __quantile_h__

//------------------------------
// Includes
//------------------------------
#include "distribution.h"
#include "locationm.h"
#include <algorithm>
#include <memory>

//------------------------------
// Class definition
//------------------------------
class CQuantile: public CDistribution
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
    virtual ~CQuantile();
  
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
  	CQuantile(SEXP radMisc);

	//-------------------
	// Private Variables
	//-------------------
    vector<double> vecd;
    double dAlpha;
    CLocationM mpLocM;
};

#endif // __quantile_h__



