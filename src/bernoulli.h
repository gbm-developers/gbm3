//------------------------------------------------------------------------------
//
//  File:       bernoulli.h
//
//  Description:   bernoulli distribution class used in GBM
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef __bernoulli_h__
#define __bernoulli_h__

//------------------------------
// Includes
//------------------------------

#include "distribution.h"
#include "dataset.h"
#include "buildinfo.h"
#include <memory>

//------------------------------
// Class definition
//------------------------------
class CBernoulli : public CDistribution
{

public:
	//---------------------
	// Factory Function
	//---------------------
	static CDistribution* Create(SEXP radMisc, const CDataset& data,
											const char* szIRMeasure,
											int& cGroups, int& cTrain);

	//---------------------
	// Public destructor
	//---------------------
    virtual ~CBernoulli();

    //---------------------
    // Public Functions
    //---------------------
    void ComputeWorkingResponse(const double *adF,
				double *adZ,
				const bag& afInBag,
				unsigned long nTrain);

    double Deviance(const double *adF,
                    unsigned long cLength,
                    bool isValidationSet=false);

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
    
    double BagImprovement(const double *adF,
                          const double *adFadj,
                          const bag& afInBag,
                          double dStepSize,
                          unsigned long nTrain);

private:
    //----------------------
    // Private Constructors
    //----------------------
    CBernoulli(SEXP radMisc, const CDataset& data);

    //-------------------
    // Private Variables
    //-------------------
    vector<double> vecdNum;
    vector<double> vecdDen;
    bool fCappedPred;
};

#endif // BERNOULLI_H



