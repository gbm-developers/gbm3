//------------------------------
//
// File: distributionFactory.h
//
// Description: factory class for distribution.h classes.
//   The gbmDataContainer class is responsible for this.
//
//------------------------------

#ifndef DISTRIBUTIONFACTORY_H
#define DISTRIBUTIONFACTORY_H

//------------------------------
// Includes
//------------------------------
#include "config_structs.h"
#include "distribution.h"
#include "dataset.h"
#include <map>
#include <memory>
#include <Rcpp.h>

//------------------------------
// Class definition
//------------------------------
class DistributionFactory
{
private:
	//----------------------
	// Private Constructors
	//----------------------
	DistributionFactory(const DistributionFactory&);

	//----------------------
	// Private Operators
	//----------------------
	DistributionFactory& operator= (const DistributionFactory&);

	//-------------------
	// Private Variables
	//-------------------
	typedef CDistribution* (*distCreate)(DataDistParams&);
	std::map<std::string, distCreate> factorymap_;

public:
	//---------------------
    // Public Constructor
	//---------------------
	DistributionFactory();

	//---------------------
	// Public destructor
	//---------------------
	~DistributionFactory(){};

	//---------------------
	// Public Functions
	//---------------------
	void RegisterDist(const std::string& kDistFamily, distCreate ptr_to_dist_createfunc);
	CDistribution* CreateDist(DataDistParams& distparams);
};

#endif // DISTRIBUTIONFACTORY_H
