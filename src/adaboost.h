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
	static CDistribution* Create(DataDistParams& distParams);

	//---------------------
	// Public destructor
	//---------------------
    virtual ~CAdaBoost();

    //---------------------
    // Public Functions
    //---------------------
    void ComputeWorkingResponse(const CDataset& data, const double *adF,
							double *adZ);

    double InitF(const CDataset& data);

    void FitBestConstant(const CDataset& data, const double *adF,
			 unsigned long cTermNodes, double* adZ, CTreeComps& treeComps);
    
    double Deviance(const CDataset& data, const double *adF,
				bool isValidationSet=false);

    double BagImprovement(const CDataset& data, const double *adF,
			  const double shrinkage, const double* adFadj);

private:
    //----------------------
    // Private Constructors
    //----------------------
    CAdaBoost();

	//-------------------
	// Private Variables
	//-------------------
   vector<double> vecdNum;
   vector<double> vecdDen;
};

#endif // ADABOOST_H


