//------------------------------------------------------------------------------
//
//  File:       tdist.h
//
//  Description:   Distribution object to implement t-distribution
//
//  History:    04/04/2008   Created
//
//------------------------------------------------------------------------------

#ifndef __tdistGBM_h__
#define __tdistGBM_h__

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
class CTDist : public CDistribution
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
	virtual ~CTDist();
  
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
			 const std::vector<unsigned long> &aiNodeAssign,
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
    CTDist(SEXP radMisc);

	//-------------------
	// Private Variables
	//-------------------
    double mdNu;
    CLocationM mpLocM;
};

#endif // __tdistGBM_h__



