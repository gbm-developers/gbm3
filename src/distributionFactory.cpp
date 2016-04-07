//-----------------------------------
//
// File: distributionFactory.cpp
//
// Description: factory class that creates the relevant distributions
//   given data and the family string.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "distributionFactory.h"
#include "adaboost.h"
#include "bernoulli.h"
#include "coxph.h"
#include "gamma.h"
#include "gaussian.h"
#include "huberized.h"
#include "laplace.h"
#include "pairwise.h"
#include "poisson.h"
#include "quantile.h"
#include "tdist.h"
#include "tweedie.h"
#include <stdexcept>

//----------------------------------------
// Function Members - Public
//----------------------------------------

//-----------------------------------
// Function: DistributionFactory
//
// Returns: none
//
// Description: initializes the factory and
//   registers the distributions to the map.
//
// Parameters: none
//
//-----------------------------------
DistributionFactory::DistributionFactory()
{
	RegisterDist("adaboost", &CAdaBoost::Create);
	RegisterDist("bernoulli", &CBernoulli::Create);
	RegisterDist("coxph", &CCoxPH::Create);
	RegisterDist("gamma", &CGamma::Create);
	RegisterDist("gaussian", &CGaussian::Create);
	RegisterDist("huberized", &CHuberized::Create);
	RegisterDist("laplace", &CLaplace::Create);
	RegisterDist("pairwise", &CPairwise::Create);
	RegisterDist("poisson", &CPoisson::Create);
	RegisterDist("quantile", &CQuantile::Create);
	RegisterDist("tdist", &CTDist::Create);
	RegisterDist("tweedie", &CTweedie::Create);
}

//----------------------------------------
// Function: RegisterDist
//
// Returns: void
//
// Description: registers family string and distribution create
//   method to factory map for later use.
//
// Parameters:
//  DistFamily - string identifying the distribution to create, IN, const string.
//  PtrDistCreateFn - function pointer to distribution create method, IN,
//                    function pointer.
//
//------------------------------------------
void DistributionFactory::RegisterDist(const string& DistFamily, distCreate PtrDistCreateFn)
{
	FactoryMap[DistFamily] = PtrDistCreateFn;
}

//------------------------------------------
// Function: CreateDist
//
// Returns: auto_ptr<CDistribution>
//
// Description: given a string looks up and creates the appropriate
//              distribution.
//
// Parameters:
//   DistFamily - string identifying the distribution to create, IN, const string.
//   radMisc - adMisc object passed in from R interface, IN, SEXP.
//   szIRMeasure - IR Measure used by the Pairwise distribution, IN, const char *.
//   cGroups  - number of groups (Pairwise) in data, IN, int&.
//   cTrain  - number of data points in training set, IN, int&.
//
//---------------------------------------------

CDistribution* DistributionFactory::CreateDist(const string& DistFamily,
															 SEXP radMisc,
															 const char* szIRMeasure,
															 int& cTrain)
{
  std::map<std::string, distCreate>::iterator it = FactoryMap.find(DistFamily);
	if( it != FactoryMap.end() )
	{
		return it -> second(radMisc, szIRMeasure, cTrain);
	}
	else
	{
		throw GBM::invalid_argument( "Error: Family string provided not recognised - distribution can't be initialized.");
	}

}
