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

CDistribution* gbm_setup
(
    double *adY,
    double *adOffset,
    double *adX,
    int *aiXOrder,
    double *adWeight,
    double *adMisc,
    int cRows,
    int cCols,
    int *acVarClasses,
    int *alMonotoneVar,
    const std::string& family,
    int cTrees,
    int cDepth,
    int cMinObsInNode,
    int cNumClasses,
    double dShrinkage,
    double dBagFraction,
    int cTrain,
    int cFeatures,
    CDataset *pData,
    int& cGroups
)
{
  CDistribution* pDist = 0;
  cGroups = -1;
  
  pData->SetData(adX,aiXOrder,adY,adOffset,adWeight,adMisc,
		 cRows,cCols,acVarClasses,alMonotoneVar);
  
    // set the distribution
  if (family == "gamma") {
      pDist = new CGamma();
  } 
  else if (family == "tweedie") {
      pDist = new CTweedie(adMisc[0]);
  }
  else if (family == "bernoulli") 
    {
      pDist = new CBernoulli();
    }
  else if (family == "gaussian") 
    {
      pDist = new CGaussian();
    }
  else if (family == "poisson")
    {
      pDist = new CPoisson();
    }
  else if (family == "adaboost")
    {
      pDist = new CAdaBoost();
    }
  else if (family == "coxph")
    {
      pDist = new CCoxPH();
    }
  else if (family == "laplace")
    {
      pDist = new CLaplace();
    }
  else if (family == "quantile")
    {
      pDist = new CQuantile(adMisc[0]);
    }
  else if (family == "tdist")
    {
      pDist = new CTDist(adMisc[0]);
    }
  else if (family == "multinomial")
    {
      pDist = new CMultinomial(cNumClasses, cRows);
    }
  else if (family == "huberized")
    {
      pDist = new CHuberized();
    }
  else if (family == "pairwise_conc")
    {
      pDist = new CPairwise("conc");
    }
  else if (family == "pairwise_ndcg")
    {
      pDist = new CPairwise("ndcg");
    }
  else if (family == "pairwise_map")
    {
      pDist = new CPairwise("map");
    }
  else if (family == "pairwise_mrr")
    {
      pDist = new CPairwise("mrr");
    }
  else
    {
      throw GBM::invalid_argument();
    }

  if (0==family.compare(0, 8, "pairwise")) 
    {
      cGroups = num_groups(adMisc, cTrain);
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




