//------------------------------
//
// File: distributionFactory.h
//
// Description: factory class for distribution.h classes.
//   Implemented in a singleton pattern fashion.
//
//------------------------------

#ifndef __distributionFactory_h__
#define __distributionFactory_h__

//------------------------------
// Includes
//------------------------------
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
	DistributionFactory();
	DistributionFactory(const DistributionFactory&);

	//----------------------
	// Private Operators
	//----------------------
	DistributionFactory& operator= (const DistributionFactory&);

	//-------------------
	// Private Variables
	//-------------------
	typedef CDistribution* (*distCreate)(SEXP, const CDataset&, const char*, int&, int&);
	std::map<std::string, distCreate> FactoryMap;

public:
	//---------------------
	// Public destructor
	//---------------------
	~DistributionFactory(){};

	//---------------------
	// Public Functions
	//---------------------
	static DistributionFactory* Get()
	{
		static DistributionFactory Instance;
		return &Instance;
	}
	void RegisterDist(const std::string& DistFamily, distCreate PtrDistCreateFn);
	CDistribution* CreateDist(const std::string& DistFamily,
					SEXP radMisc, const CDataset& data,
					const char* szIRMeasure, int& cGroups, int& cTrain);
};

#endif // __distributionFactory_h__
