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
#include <vector>
#include <memory>
#include "buildinfo.h"
#include "distribution.h"
#include "tree.h"
#include "dataset.h"
#include "node_factory.h"

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
    void Initialize(CDistribution *pDist,
		    double dLambda,
		    unsigned long nTrain,
		    unsigned long cFeatures,
		    double dBagFraction,
		    unsigned long cLeaves,
		    unsigned long cMinObsInNode,
		    int cGroups);

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

    bool IsPairwise() const { return (cGroups >= 0); }

private:
	//-------------------
	// Private Variables
	//-------------------
    CDistribution *pDist;       // the distribution - this contains the data
    bool fInitialized;          // indicates whether the GBM has been initialized
    std::auto_ptr<CNodeFactory> pNodeFactory;

    // these objects are for the tree growing
    // allocate them once here for all trees to use
    bag afInBag;
    std::vector<unsigned long> aiNodeAssign;
    std::vector<CNodeSearch> aNodeSearch;
    std::auto_ptr<CCARTTree> ptreeTemp;
    VEC_P_NODETERMINAL vecpTermNodes;
    std::vector<double> adZ;
    std::vector<double> adFadj;

    double dLambda;
    unsigned long cTrain;
    unsigned long cValid;
    unsigned long cFeatures;
    unsigned long cTotalInBag;
    double dBagFraction;
    unsigned long cDepth;
    unsigned long cMinObsInNode;
    int  cGroups;
};

#endif //  __gbm_enginegbm_h__



