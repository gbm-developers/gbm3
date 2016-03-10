//------------------------------------------------------------------------------
//
//  File:       coxph.h
//
//  Description:   Cox proportional hazard object
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef __coxph_h__
#define __coxph_h__

//------------------------------
// Includes
//------------------------------

#include "distribution.h"
#include "matrix.h"
#include <memory>

//------------------------------
// Class definition
//------------------------------
class CCoxPH : public CDistribution
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
    virtual ~CCoxPH();

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
    CCoxPH(SEXP radMisc);

    //-------------------
    // Private Variables
    //-------------------
    vector<double> vecdP;
    vector<double> vecdRiskTot;
    vector<double> vecdG;
    vector<unsigned long> veciK2Node;
    vector<unsigned long> veciNode2K;

    matrix<double> matH;
    matrix<double> matHinv;
    const double* adDelta;
};

#endif // __coxph_h__



