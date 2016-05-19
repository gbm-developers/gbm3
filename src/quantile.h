//------------------------------------------------------------------------------
//  File:       quantile.h
//
//  Contents:   quantile regression for GBM.
//
//  History:    10/8/2006   Created by Brian Kriegler (bk@stat.ucla.edu)
//              6/11/2007   gregr merged with official gbm
//
//------------------------------------------------------------------------------

#ifndef QUANTILE_H
#define QUANTILE_H

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
class CQuantile: public CDistribution
{

public:
	//---------------------
	// Factory Function
	//---------------------
	static CDistribution* Create(DataDistParams& distParams);

	//---------------------
	// Public destructor
	//---------------------
    virtual ~CQuantile();
  
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
  	CQuantile(double alpha);

	//-------------------
	// Private Variables
	//-------------------
    vector<double> vecd;
    double dAlpha;
    CLocationM mpLocM;
};

#endif // QUANTILE_H



