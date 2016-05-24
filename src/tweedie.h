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
	static CDistribution* Create(DataDistParams& distParams);

	//---------------------
	// Public destructor
	//---------------------
    virtual ~CTweedie();

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
			  const double shrinkage,
                          const double* adFadj);
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

