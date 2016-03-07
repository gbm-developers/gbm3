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
#include "gbm_setup.h"
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
    CGBM();

	//---------------------
	// Public destructor
	//---------------------
    ~CGBM();

	//---------------------
	// Public Functions
	//---------------------
    void Initialize();
    void SetDataAndDistribution(const CDataset& data, SEXP radMisc, const std::string& family,
    		const int cTrain, int& cGroups);
    void SetTreeContainer(double dLambda,
    	    unsigned long cTrain,
    	    unsigned long cFeatures,
    	    double dBagFraction,
    	    unsigned long cDepth,
    	    unsigned long cMinObsInNode,
    	    int cGroups);

    void iterate(double *adF,
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

    void printDist(){
    	std::cout << pDist->data_ptr()->nrow() << endl;
    }

private:
	//-------------------
	// Private Variables
	//-------------------
    CDistribution* pDist;       // the distribution - this contains the data
    CNodeFactory* pNodeFactory;
    bool fInitialized;          // indicates whether the GBM has been initialized

    // these objects are for the tree growing
    // allocate them once here for all trees to use
    CTreeComps* pTreeComp;
};

#endif //  __gbm_enginegbm_h__



