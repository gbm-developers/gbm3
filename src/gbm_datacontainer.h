//------------------------------------------------------------------------------
//
 //  File:       gbmDataContainer.h
//
//  Description:   Header file for a class that stores the data and dist for gbm.
//
//------------------------------------------------------------------------------


#ifndef GBMDATACONTAINER_H
#define GBMDATACONTAINER_H

//------------------------------
// Includes
//------------------------------
#include "buildinfo.h"
#include "config_structs.h"
#include "dataset.h"
#include "distribution.h"
#include "distribution_factory.h"
#include "gbm_treecomponents.h"
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
    CGBMDataContainer(DataDistParams& datadist_config);

	//---------------------
	// Public destructor
	//---------------------
    ~CGBMDataContainer();

    //---------------------
	// Public Functions
	//---------------------
    void Initialize();
    double InitialFunctionEstimate();
    void ComputeResiduals(const double* kFuncEstimate, double* residuals);
    void ComputeBestTermNodePreds(const double* kFuncEstimate, double* residuals, CTreeComps& treecomps);
    double ComputeDeviance(const double *kFuncEstimate, bool is_validationset=false);
    double ComputeBagImprovement(const double* kFuncEstimate, const double kShrinkage, const double* kDeltaEstimate);
    void BagData();
    
    CDistribution* get_dist(){ return distptr_; }
    const CDataset& get_data() const { return data_; }
    CDataset& get_data() { return data_; }

private:
	//-------------------
	// Private Variables
	//-------------------
    CDataset data_;
    CDistribution* distptr_;
    DistributionFactory* distfactory_;
};

#endif // GBMDATACONTAINER_H
