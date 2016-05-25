//------------------------------------------------------------------------------
//
//  File:       tweedie.h
//
//  Description: tweedie distribution for .
//
//------------------------------------------------------------------------------

#ifndef TWEEDIE_H
#define TWEEDIE_H

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
	static CDistribution* Create(DataDistParams& distparams);

	//---------------------
	// Public destructor
	//---------------------
    virtual ~CTweedie();

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
    				const double* kFuncEstimate);
    
    double BagImprovement(const CDataset& kData,
			  const double* kFuncEstimate,
			  const double kShrinkage,
                          const double* kDeltaEstimates);
private:
    //----------------------
    // Private Constructors
    //----------------------
    CTweedie(double power);

	//-------------------
	// Private Variables
	//-------------------
    double power_;
};

#endif // TWEEDIE_H

