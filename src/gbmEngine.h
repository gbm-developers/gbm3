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
#include "gbmDataContainer.h"
#include "node_factory.h"
#include "gbmTreeComps.h"
#include <memory>
#include <Rcpp.h>
#include <vector>

//------------------------------
// Forward declaration
//------------------------------
class CGBMDataContainer;

//------------------------------
// Class definition
//------------------------------
class CGBM
{
public:
	//----------------------
	// Public Constructors
	//----------------------
    CGBM();

	//---------------------
	// Public destructor
	//---------------------
    ~CGBM();

	//---------------------
	// Public Functions
	//---------------------
    void Initialize();
    void SetDataAndDistribution(SEXP radY, SEXP radOffset, SEXP radX, SEXP raiXOrder,
            SEXP radWeight, SEXP racVarClasses,
            SEXP ralMonotoneVar, SEXP radMisc, const std::string& family,
    		const int cTrain, int& cGroups);
    void SetTreeContainer(double dLambda,
    	    unsigned long cTrain,
    	    unsigned long cFeatures,
    	    double dBagFraction,
    	    unsigned long cDepth,
    	    unsigned long cMinObsInNode,
    	    int cGroups);

    void Iterate(double *adF,
		 double &dTrainError,
		 double &dValidError,
		 double &dOOBagImprove,
		 int &cNodes);
    
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

    void InitF(double &dInitF, unsigned long cLength);
    void UpdateParams(const double *adF,
            			      unsigned long cLength);
    bool IsPairwise() const { return (pTreeComp->GetNoGroups() >= 0); }

private:
	//-------------------
	// Private Variables
	//-------------------
    CGBMDataContainer* pDataCont;
    CNodeFactory* pNodeFactory;
    CTreeComps* pTreeComp;
    bool fInitialized;          // indicates whether the GBM has been initialized
    bool hasDataAndDist, hasTreeContainer;
};

#endif //  __gbmEnginegbm_h__



