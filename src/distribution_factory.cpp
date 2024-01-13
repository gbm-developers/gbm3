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
#include "distribution_factory.h"
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
DistributionFactory::DistributionFactory() {
  RegisterDist("adaboost", &CAdaBoost::Create);
  RegisterDist("bernoulli", &CBernoulli::Create);
  RegisterDist("coxph", &CCoxPH::Create);
  RegisterDist("gamma", &CGamma::Create);
  RegisterDist("gaussian", &CGaussian::Create);
  RegisterDist("huberized", &CHuberized::Create);
  RegisterDist("laplace", &CLaplace::Create);
  RegisterDist("pairwise_conc", &CPairwise::Create);
  RegisterDist("pairwise_ndcg", &CPairwise::Create);
  RegisterDist("pairwise_map", &CPairwise::Create);
  RegisterDist("pairwise_mrr", &CPairwise::Create);
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
//  DistFamily - string identifying the distribution to create, IN, const
//  string.
//  PtrDistCreateFn - function pointer to distribution create method, IN,
//                    function pointer.
//
//------------------------------------------
void DistributionFactory::RegisterDist(const std::string& kDistFamily,
                                       DistCreate ptr_to_dist_createfunc) {
  factorymap_.insert(
    std::pair<std::string, DistCreate>(kDistFamily, ptr_to_dist_createfunc));
}

//------------------------------------------
// Function: CreateDist
//
// Returns: ptr CDistribution
//
// Description: given a string looks up and creates the appropriate
//              distribution.
//
// Parameters:
//   distParams - configuration struct containing parameters for distribution
//   set-up.
//---------------------------------------------

CDistribution* DistributionFactory::CreateDist(DataDistParams& distparams) {
  std::multimap<std::string, DistCreate>::iterator it =
      factorymap_.find(distparams.family);
  if (it != factorymap_.end()) {
    return it->second(distparams);
  } else {
    throw gbm_exception::InvalidArgument(
        "Error: Family string provided not recognised - distribution can't be "
        "initialized.");
  }
}
