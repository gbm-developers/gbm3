//------------------------------------------------------------------------------
//
//  File:       gbm_engine.h
//
//  Description:   Header file for Gradient Boosting Engine.
//
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef __gbmEnginegbm_h__
#define __gbmEnginegbm_h__

//------------------------------
// Includes
//------------------------------
#include "buildinfo.h"
#include "configStructs.h"
#include "gbmDataContainer.h"
#include "gbmTreeComps.h"
#include <memory>
#include <Rcpp.h>
#include <vector>

//------------------------------
// Class definition
//------------------------------
class CGBM
{
public:
	//----------------------
	// Public Constructors
	//----------------------
    CGBM(configStructs GBMParams);

	//---------------------
	// Public destructor
	//---------------------
    ~CGBM();

	//---------------------
	// Public Functions
	//---------------------
    void FitLearner(double *adF,
		 double &dTrainError,
		 double &dValidError,
		 double &dOOBagImprove);

    void GBMTransferTreeToRList(int *aiSplitVar,
			     double *adSplitPoint,
			     int *aiLeftNode,
			     int *aiRightNode,
			     int *aiMissingNode,
			     double *adErrorReduction,
			     double *adWeight,
			     double *adPred,
			     VEC_VEC_CATEGORIES &vecSplitCodes,
			     int cCatSplitsOld);

    const long SizeOfFittedTree() const;
    double InitF();

private:
	//-------------------
	// Private Variables
	//-------------------
    CGBMDataContainer* pDataCont;
	CTreeComps* pTreeComp;
    bool fInitialized;          // indicates whether the GBM has been initialized

    // Residuals and adjustments to function estimate
    std::vector<double> adZ;

};

#endif //  __gbmEnginegbm_h__



