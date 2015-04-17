//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//
//  File:       gbm_engine.h
//
//  License:    GNU GPL (version 2 or later)
//
//  Contents:   Generalized boosted model engine
//
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef GBM_ENGINGBM_H
#define GBM_ENGINGBM_H

#include <vector>
#include <memory>
#include "buildinfo.h"
#include "distribution.h"
#include "tree.h"
#include "dataset.h"
#include "node_factory.h"

using namespace std;

class CGBM
{

public:

    CGBM();
    ~CGBM();
    void Initialize(CDataset &pData,
		    CDistribution *pDist,
		    double dLambda,
		    unsigned long nTrain,
		    unsigned long cFeatures,
		    double dBagFraction,
		    unsigned long cLeaves,
		    unsigned long cMinObsInNode,
		    unsigned long cNumClasses,
		    int cGroups);

    void iterate(double *adF,
		 double &dTrainError,
		 double &dValidError,
		 double &dOOBagImprove,
		 int &cNodes,
		 int cNumClasses,
		 int cClassIdx);
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

    CDataset *pData;            // the data
    CDistribution *pDist;       // the distribution
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

#endif // GBM_ENGINGBM_H



