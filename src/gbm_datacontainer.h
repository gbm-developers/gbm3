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
    CGBMDataContainer(DataDistParams& dataDistConfig);

	//---------------------
	// Public destructor
	//---------------------
    ~CGBMDataContainer();

    //---------------------
	// Public Functions
	//---------------------
    void Initialize();
    double InitialFunctionEstimate();
    void ComputeResiduals(const double* adF, double* adZ);
    void ComputeBestTermNodePreds(const double* adF, double* adZ, CTreeComps& pTreeComp);
    double ComputeDeviance(const double *adF, bool isValidationSet=false);
    double ComputeBagImprovement(const double* adF, const double shrinkage, const double* adFadj);
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
