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

#ifndef __poisson_h__
#define __poisson_h__

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
	static CDistribution* Create(SEXP radMisc,
										const char* szIRMeasure,
										int& cGroups, int& cTrain);

	//---------------------
	// Public destructor
	//---------------------
    virtual ~CPoisson();

    //---------------------
    // Public Functions
    //---------------------
    void ComputeWorkingResponse(const CDataset* pData,
    		const double *adF,
				double *adZ,
				const bag& afInBag,
				unsigned long nTrain);

    double Deviance(const CDataset* pData,
    				const double *adF,
                    unsigned long cLength,
                    bool isValidationSet=false);

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
    CPoisson(SEXP radMisc);


	//-------------------
	// Private Variables
	//-------------------
    vector<double> vecdNum;
    vector<double> vecdDen;
    vector<double> vecdMax;
    vector<double> vecdMin;
};

#endif // __poisson_h__



