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

#ifndef __gbm_enginegbm_h__
#define __gbm_enginegbm_h__

//------------------------------
// Includes
//------------------------------
#include "buildinfo.h"
#include "distribution.h"
#include "tree.h"
#include "dataset.h"
#include "node_factory.h"
#include "gbmTreeComps.h"
#include <vector>
#include <memory>

//------------------------------
// Class definition
//------------------------------
class CGBM
{
public:
	//----------------------
	// Public Constructors
	//----------------------
    CGBM(CDistribution* DistPtr, double dLambda,
    	    unsigned long cTrain,
    	    unsigned long cFeatures,
    	    double dBagFraction,
    	    unsigned long cDepth,
    	    unsigned long cMinObsInNode,
    	    int cGroups);


	//---------------------
	// Public destructor
	//---------------------
    ~CGBM();

	//---------------------
	// Public Functions
	//---------------------
    void Initialize();

    void iterate(double *adF,
		 double &dTrainError,
		 double &dValidError,
		 double &dOOBagImprove,
		 int &cNodes);
    
    void TransferTreeToRList(int *aiSplitVar,
			     double *adSplitPoint,
			     int *aiLeftNode,
			     int *aiRightNode,
			     int *aiMissingNode,
			     double *adErrorReduction,
			     double *adWeight,
			     double *adPred,
			     VEC_VEC_CATEGORIES &vecSplitCodes,
			     int cCatSplitsOld);

    bool IsPairwise() const { return (pTreeComp->GetNoGroups() >= 0); }

private:
	//-------------------
	// Private Variables
	//-------------------
    CDistribution* pDist;       // the distribution - this contains the data
    bool fInitialized;          // indicates whether the GBM has been initialized
    std::auto_ptr<CNodeFactory> pNodeFactory;

    // these objects are for the tree growing
    // allocate them once here for all trees to use
    CTreeComps* pTreeComp;
};

#endif //  __gbm_enginegbm_h__



