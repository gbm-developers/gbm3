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
										const char* szIRMeasure,
										int& cGroups, int& cTrain);

	//---------------------
	// Public destructor
	//---------------------
    virtual ~CTweedie();

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
    CTweedie(SEXP radMisc);

	//-------------------
	// Private Variables
	//-------------------
    vector<double> vecdNum;
    vector<double> vecdDen;
    vector<double> vecdMax;
    vector<double> vecdMin;
    double dPower;
};

#endif // _tweedie_h__

