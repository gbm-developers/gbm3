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
										const char* szIRMeasure,
										int& cGroups, int& cTrain);

	//---------------------
	// Public destructor
	//---------------------
	virtual ~CLaplace();

	//---------------------
	// Public Functions
	//---------------------
	void ComputeWorkingResponse(const CDataset* pData,
				const double *adF,
				  double *adZ,
				  const bag& afInBag,
				  unsigned long nTrain);

	void InitF(const CDataset* pData,
			double &dInitF,
			unsigned long cLength);
  
	void FitBestConstant(const CDataset* pData,
				   const double *adF,
		       	   double *adZ,
		       	   const std::vector<unsigned long>& aiNodeAssign,
		       	   unsigned long nTrain,
		       	   VEC_P_NODETERMINAL vecpTermNodes,
		       	   unsigned long cTermNodes,
		       	   unsigned long cMinObsInNode,
		       	   const bag& afInBag,
		       	   const double *adFadj);
  
	double Deviance(const CDataset* pData,
					const double *adF,
                    unsigned long cLength,
                    bool isValidationSet=false);
  
	double BagImprovement(const CDataset* pData,
				const double *adF,
				const double *adFadj,
				const bag& afInBag,
				double dStepSize,
				unsigned long nTrain);

private:
	//----------------------
	// Private Constructors
	//----------------------
	CLaplace(SEXP radMisc);

	//-------------------
	// Private Variables
	//-------------------
	vector<double> vecd;
	vector<double>::iterator itMedian;
	CLocationM mpLocM;
};

#endif // __laplaceGBM_h__



