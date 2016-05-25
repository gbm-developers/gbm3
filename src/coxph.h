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

#ifndef COXPH_H
#define COXPH_H

//------------------------------
// Includes
//------------------------------

#include "distribution.h"
#include <memory>

//------------------------------
// Class Forwards and Enums
//------------------------------
class GenericCoxState;

//------------------------------
// Class definition
//------------------------------
class CCoxPH : public CDistribution
{

public:
	//---------------------
	// Factory Function
	//---------------------
	static CDistribution* Create(DataDistParams& distparams);

	//---------------------
	// Public destructor
	//---------------------
    virtual ~CCoxPH();

    //---------------------
    // Public Functions
    //---------------------
    void ComputeWorkingResponse(const CDataset& kData,
    		const double* kFuncEstimate,
				double* residuals);

    double InitF(const CDataset& kData);
    
    void FitBestConstant(const CDataset& kData,
    		const double* kFuncEstimate,
			 unsigned long num_terminalnodes,
			 double* residuals,
			 CCARTTree& tree);
    
    double Deviance(const CDataset& kData,
    				const double *kFuncEstimate);

    double BagImprovement(const CDataset& kData,
			  const double* kFuncEstimate,
			  const double kShrinkage, const double* kDeltaEstimate);
    void ShiftDistPtrs(unsigned long shift)
    {
    	status_ = shift_ptr(status_, shift);
		sortedendtimes_ = shift_ptr(sortedendtimes_, shift);
		sortedstarttimes_ = shift_ptr(sortedstarttimes_, shift);
		strata = shift_ptr(strata, shift);
    }

    // Getters for the internal variables
    double* StatusVec();
    const double* StatusVec() const;

    int* EndTimeIndices();
    const int* EndTimeIndices() const;

    int* StartTimeIndices();
    const int* StartTimeIndices() const;

    int* StrataVec();
    const int* StrataVec() const;

    int TieApproxMethod() const;

    double PriorCoeffVar() const;

private:
    //----------------------
    // Private Constructors
    //----------------------
    CCoxPH(double* stats, int* sorted_end, int* sorted_start, int* strats,
	   bool is_startstop, int tiesmethod, double priorcoeff);

    //----------------------
    // Private Functions
    //----------------------
    double LogLikelihood(const int n, const CDataset& kData,
			 const double* eta, double* resid);

    double LogLikelihoodTiedTimes(const int n, const CDataset& kData,
				  const double* eta, double* resid);

    //-------------------
    // Private Variables
    //-------------------
    const bool kStartStopCase_;
    int* sortedendtimes_;
    int* sortedstarttimes_;
    int* strata;
    const double kPriorCoeffVariation_;
    double* status_;
    int tiedtimesmethod_;
    GenericCoxState* coxstate_methods_;

    
  	
  	
};

#endif // COXPH_H



