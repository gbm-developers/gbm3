//------------------------------------------------------------------------------
//
//  File:       gbmDataContainer.h
//
//  Description:   Header file for a class that stores the data and dist for gbm.
//
//------------------------------------------------------------------------------


#ifndef __gbmDataContainer_h__
#define __gbmDataContainer_h__

//------------------------------
// Includes
//------------------------------
#include "buildinfo.h"
#include "dataset.h"
#include "distribution.h"
#include "distributionFactory.h"
#include "gbmTreeComps.h"
#include <Rcpp.h>
#include <vector>
#include <memory>

//------------------------------
// Class definition
//------------------------------
class CGBMDataContainer
{
public:
	//----------------------
	// Public Constructors
	//----------------------
    CGBMDataContainer(SEXP radY, SEXP radOffset, SEXP radX, SEXP raiXOrder,
            SEXP radWeight, SEXP racVarClasses,
            SEXP ralMonotoneVar, SEXP radMisc, const std::string& family, int cTrain, int& cGroups);

	//---------------------
	// Public destructor
	//---------------------
    ~CGBMDataContainer();

    //---------------------
	// Public Functions
	//---------------------
    void Initialize();
    void InitializeFunctionEstimate(double &dInitF, unsigned long cLength);
    void ComputeResiduals(const double* adF, CTreeComps* pTreeComp);
    void ComputeBestTermNodePreds(const double* adF, CTreeComps* pTreeComp, int& cNodes);
    double ComputeDeviance(const double *adF, CTreeComps* pTreeComp,  bool isValidationSet=false);
    double ComputeBagImprovement(const double* adF, CTreeComps* pTreeComp);
    CDistribution* getDist();
    const CDataset* getData();


private:
	//-------------------
	// Private Variables
	//-------------------
    const CDataset data;
    CDistribution* pDist;
    DistributionFactory* DistFactory; // currently a singleton - does not need to be now - TODO:remove property.

};

#endif //  __gbmDataContainer_h__
