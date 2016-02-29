//------------------------------------------------------------------------------
//
//  File:       gaussian.h
//
//  Contents:   gaussian distribution for GBM.
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//------------------------------------------------------------------------------

#ifndef __gaussian_h__
#define __gaussian_h__

//------------------------------
// Includes
//------------------------------
#include "distribution.h"
#include "dataset.h"
#include <memory>

//------------------------------
// Class definition
//------------------------------
class CGaussian : public CDistribution
{

public:
	//---------------------
	// Factory Function
	//---------------------
	static std::auto_ptr<CDistribution> Create(SEXP radMisc, const CDataset& data,
											const char* szIRMeasure,
											int& cGroups, int& cTrain);

	//---------------------
	// Public destructor
	//---------------------
    virtual ~CGaussian();

    //---------------------
    // Public Functions
    //---------------------
    void ComputeWorkingResponse(const double *adF,
				double *adZ,
				const bag& afInBag,
				unsigned long nTrain);

    void InitF(double &dInitF,
	       unsigned long cLength);
    
    void FitBestConstant(const double *adF,
					 double *adZ,
					 const std::vector<unsigned long>& aiNodeAssign,
					 unsigned long nTrain,
					 VEC_P_NODETERMINAL vecpTermNodes,
					 unsigned long cTermNodes,
					 unsigned long cMinObsInNode,
					 const bag& afInBag,
					 const double *adFadj);

    double Deviance(const double *adF,
                    unsigned long cLength,
                    bool isValidationSet=false);

    double BagImprovement(const double *adF,
                          const double *adFadj,
                          const bag& afInBag,
                          double dStepSize,
                          unsigned long nTrain);
private:
    //----------------------
    // Private Constructors
    //----------------------
    CGaussian(SEXP radMisc, const CDataset& data);
};

#endif // __gaussian_h__



