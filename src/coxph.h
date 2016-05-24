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
	static CDistribution* Create(DataDistParams& distParams);

	//---------------------
	// Public destructor
	//---------------------
    virtual ~CCoxPH();

    //---------------------
    // Public Functions
    //---------------------
    void ComputeWorkingResponse(const CDataset& data,
    		const double *adF,
				double *adZ);

    double InitF(const CDataset& data);
    
    void FitBestConstant(const CDataset& data,
    		const double *adF,
			 unsigned long cTermNodes,
			 double* adZ,
			 CTreeComps& treeComps);
    
    double Deviance(const CDataset& data,
    				const double *adF,
                    bool isValidationSet=false);

    double BagImprovement(const CDataset& data,
			  const double *adF,
			  const double shrinkage, const double* adFadj);

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
    CCoxPH(double* stats, int* sortedEnd, int* sortedSt, int* strats,
	   bool isStartStop, int tiedMethod, double priorCoeff);

    //----------------------
    // Private Functions
    //----------------------
    double LogLikelihood(const int n, const CDataset& data,
			 const double* eta, double* resid);

    double LogLikelihoodTiedTimes(const int n, const CDataset& data,
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



