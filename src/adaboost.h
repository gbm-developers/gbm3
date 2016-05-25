//------------------------------------------------------------------------------
//
//  File:       adaboost.h
//
//  Description: distribution used in adaboost fitting.
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef ADABOOST_H
#define ADABOOST_H

//------------------------------
// Includes
//------------------------------
#include "distribution.h"
#include <memory>

//------------------------------
// Class definition
//------------------------------
class CAdaBoost : public CDistribution
{

public:
	//---------------------
	// Factory Function
	//---------------------
	static CDistribution* Create(DataDistParams& distparams);

	//---------------------
	// Public destructor
	//---------------------
    virtual ~CAdaBoost();

    //---------------------
    // Public Functions
    //---------------------
    void ComputeWorkingResponse(const CDataset& kData, const double* kFuncEstimate,
							double* residuals);

    double InitF(const CDataset& kData);

    void FitBestConstant(const CDataset& kData, const double* kFuncEstimate,
			 unsigned long numterminal_nodes, double* residuals, CCARTTree& tree);
    
    double Deviance(const CDataset& kData, const double* kFuncEstimate);

    double BagImprovement(const CDataset& kData, const double* kFuncEstimate,
			  const double shrinkage, const double* kDeltaEstimate);

private:
    //----------------------
    // Private Constructors
    //----------------------
    CAdaBoost();

	//-------------------
	// Private Variables
	//-------------------
   vector<double> numerator_bestconstant_;
   vector<double> denominator_bestconstant_;
};

#endif // ADABOOST_H


