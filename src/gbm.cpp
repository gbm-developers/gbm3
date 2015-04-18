//------------------------------------------------------------------------------
//
//  GBM by Greg Ridgeway  Copyright (C) 2003
//  File:       gbm.cpp
//
//------------------------------------------------------------------------------

#include <algorithm>
#include <vector>
#include "gbm.h"

// Count the number of distinct groups in the input data
int num_groups(const double* adMisc, int cTrain)
{
    if (cTrain <= 0)
    {
        return 0;
    }
    double dLastGroup = adMisc[0];
    int cGroups = 1;

    for(int i=1; i<cTrain; i++)
    {
        const double dGroup = adMisc[i];
        if (dGroup != dLastGroup)
        {
            dLastGroup = dGroup;
            cGroups++;
        }
    }
    return cGroups;
}

std::auto_ptr<CDistribution> gbm_setup
(
 const CDataset& data,
 const std::string& family,
 int cTrees,
 int cDepth,
 int cMinObsInNode,
 int cNumClasses,
 double dShrinkage,
 double dBagFraction,
 int cTrain,
 int cFeatures,
 int& cGroups
 )
{
  std::auto_ptr<CDistribution> pDist;
  cGroups = -1;
  
    // set the distribution
  if (family == "gamma") {
    pDist.reset(new CGamma());
  } 
  else if (family == "tweedie") {
    pDist.reset(new CTweedie(data.misc_ptr()[0]));
  }
  else if (family == "bernoulli") 
    {
      pDist.reset(new CBernoulli());
    }
  else if (family == "gaussian") 
    {
      pDist.reset(new CGaussian());
    }
  else if (family == "poisson")
    {
      pDist.reset(new CPoisson());
    }
  else if (family == "adaboost")
    {
      pDist.reset(new CAdaBoost());
    }
  else if (family == "coxph")
    {
      pDist.reset(new CCoxPH());
    }
  else if (family == "laplace")
    {
      pDist.reset(new CLaplace());
    }
  else if (family == "quantile")
    {
      pDist.reset(new CQuantile(data.misc_ptr()[0]));
    }
  else if (family == "tdist")
    {
      pDist.reset(new CTDist(data.misc_ptr()[0]));
    }
  else if (family == "multinomial")
    {
      pDist.reset(new CMultinomial(cNumClasses, data.nrow()));
    }
  else if (family == "huberized")
    {
      pDist.reset(new CHuberized());
    }
  else if (family == "pairwise_conc")
    {
      pDist.reset(new CPairwise("conc"));
    }
  else if (family == "pairwise_ndcg")
    {
      pDist.reset(new CPairwise("ndcg"));
    }
  else if (family == "pairwise_map")
    {
      pDist.reset(new CPairwise("map"));
    }
  else if (family == "pairwise_mrr")
    {
      pDist.reset(new CPairwise("mrr"));
    }
  else
    {
      throw GBM::invalid_argument();
    }

  if (0==family.compare(0, 8, "pairwise")) 
    {
      cGroups = num_groups(data.misc_ptr(), cTrain);
    }
  
  return pDist;
}


void gbm_transfer_to_R
(
 CGBM *pGBM,
 VEC_VEC_CATEGORIES &vecSplitCodes,
 int *aiSplitVar,
 double *adSplitPoint,
 int *aiLeftNode,
 int *aiRightNode,
 int *aiMissingNode,
 double *adErrorReduction,
 double *adWeight,
 double *adPred,
 int cCatSplitsOld
 )
{
    pGBM->TransferTreeToRList(aiSplitVar,
			      adSplitPoint,
			      aiLeftNode,
			      aiRightNode,
			      aiMissingNode,
			      adErrorReduction,
			      adWeight,
			      adPred,
			      vecSplitCodes,
			      cCatSplitsOld);
}




