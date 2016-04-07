//------------------------------------------------------------------------------
//
//  File:       adaboost.h
//
//  Description: distribution used in adaboost fitting.
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef __adaboost_h__
#define __adaboost_h__

//------------------------------
// Includes
//------------------------------
#include "distribution.h"
#include <memory>

//------------------------------
// Class definition
//------------------------------
class CAdaBoost : public CDistribution
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
    virtual ~CAdaBoost();

    //---------------------
    // Public Functions
    //---------------------
    void ComputeWorkingResponse(const CDataset* pData, const double *adF,
							double *adZ);

    double InitF(const CDataset* pData);

    void FitBestConstant(const CDataset* pData, const double *adF,
					 unsigned long cTermNodes, double* adZ, CTreeComps* pTreeComps);
    
    double Deviance(const CDataset* pData, const double *adF,
				bool isValidationSet=false);

    double BagImprovement(const CDataset& data, const double *adF,
    						const bag& afInBag,
					  	  const double shrinkage, const double* adFadj);

private:
    //----------------------
    // Private Constructors
    //----------------------
    CAdaBoost(SEXP radMisc);

	//-------------------
	// Private Variables
	//-------------------
   vector<double> vecdNum;
   vector<double> vecdDen;
};

#endif // __adaboost_h__


