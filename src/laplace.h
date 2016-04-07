//------------------------------------------------------------------------------
//
//  File:       laplace.h
//
//  Description:   laplace distribution for GBM.
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//------------------------------------------------------------------------------

#ifndef __laplaceGBM_h__
#define __laplaceGBM_h__

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
class CLaplace : public CDistribution
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
	virtual ~CLaplace();

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
	CLaplace(SEXP radMisc);

	//-------------------
	// Private Variables
	//-------------------
	CLocationM mpLocM;
};

#endif // __laplaceGBM_h__



