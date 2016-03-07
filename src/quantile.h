//------------------------------------------------------------------------------
//  File:       quantile.h
//
//  Contents:   quantile regression for GBM.
//
//  History:    10/8/2006   Created by Brian Kriegler (bk@stat.ucla.edu)
//              6/11/2007   gregr merged with official gbm
//
//------------------------------------------------------------------------------

#ifndef __quantile_h__
#define __quantile_h__

//------------------------------
// Includes
//------------------------------
#include "distribution.h"
#include "dataset.h"
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
	static CDistribution* Create(SEXP radMisc, const CDataset& data,
										const char* szIRMeasure,
										int& cGroups, int& cTrain);

	//---------------------
	// Public destructor
	//---------------------
    virtual ~CQuantile();
  
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
  	CQuantile(SEXP radMisc, const CDataset& data);

	//-------------------
	// Private Variables
	//-------------------
    vector<double> vecd;
    double dAlpha;
    CLocationM mpLocM;
};

#endif // __quantile_h__



