//------------------------------------------------------------------------------
//
//  GBM by Greg Ridgeway  Copyright (C) 2003
//  File:       gbm.cpp
//
//------------------------------------------------------------------------------

#include <algorithm>
#include<vector>
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

unsigned long gbm_setup
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
    const char *pszFamily,
    int cTrees,
    int cDepth,
    int cMinObsInNode,
    int cNumClasses,
    double dShrinkage,
    double dBagFraction,
    int cTrain,
	int cFeatures,
    CDataset *pData,
    PCDistribution &pDist,
    int& cGroups
)
{
    unsigned long hr = 0;
    cGroups = -1;

    hr = pData->SetData(adX,aiXOrder,adY,adOffset,adWeight,adMisc,
                        cRows,cCols,acVarClasses,alMonotoneVar);

    if(GBM_FAILED(hr))
    {
        goto Error;
    }

    // set the distribution
    if(strncmp(pszFamily,"bernoulli",2) == 0)
    {
        pDist = new CBernoulli();
    }
    else if(strncmp(pszFamily,"gaussian",2) == 0)
    {
        pDist = new CGaussian();
    }
    else if(strncmp(pszFamily,"poisson",2) == 0)
    {
        pDist = new CPoisson();
    }
    else if(strncmp(pszFamily,"adaboost",2) == 0)
    {
        pDist = new CAdaBoost();
    }
    else if(strncmp(pszFamily,"coxph",2) == 0)
    {
        pDist = new CCoxPH();
    }
    else if(strncmp(pszFamily,"laplace",2) == 0)
    {
        pDist = new CLaplace();
    }
    else if(strncmp(pszFamily,"quantile",2) == 0)
    {
        pDist = new CQuantile(adMisc[0]);
    }
    else if(strncmp(pszFamily,"tdist",2) == 0)
    {
        pDist = new CTDist(adMisc[0]);
    }
    else if(strncmp(pszFamily,"multinomial",2) == 0)
    {
        pDist = new CMultinomial(cNumClasses, cRows);
    }
    else if(strncmp(pszFamily,"huberized",2) == 0)
    {
        pDist = new CHuberized();
    }
    else if(strcmp(pszFamily,"pairwise_conc") == 0)
    {
        pDist = new CPairwise("conc");
    }
    else if(strcmp(pszFamily,"pairwise_ndcg") == 0)
    {
        pDist = new CPairwise("ndcg");
    }
    else if(strcmp(pszFamily,"pairwise_map") == 0)
    {
        pDist = new CPairwise("map");
    }
    else if(strcmp(pszFamily,"pairwise_mrr") == 0)
    {
        pDist = new CPairwise("mrr");
    }
    else
    {
        hr = GBM_INVALIDARG;
        goto Error;
    }

    if (!strncmp(pszFamily, "pairwise", strlen("pairwise")))
    {
        cGroups = num_groups(adMisc, cTrain);
    }

Cleanup:
    return hr;
Error:
    goto Cleanup;
}


GBMRESULT gbm_transfer_to_R
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
    GBMRESULT hr = GBM_OK;


    hr = pGBM->TransferTreeToRList(aiSplitVar,
                                   adSplitPoint,
                                   aiLeftNode,
                                   aiRightNode,
                                   aiMissingNode,
                                   adErrorReduction,
                                   adWeight,
                                   adPred,
                                   vecSplitCodes,
                                   cCatSplitsOld);
    if(GBM_FAILED(hr)) goto Error;

Cleanup:
    return hr;
Error:
    goto Cleanup;
}


GBMRESULT gbm_transfer_catsplits_to_R
(
    int iCatSplit,
    VEC_VEC_CATEGORIES &vecSplitCodes,
    int *aiSplitCodes
)
{
    std::copy(vecSplitCodes[iCatSplit].begin(),
	      vecSplitCodes[iCatSplit].end(),
	      aiSplitCodes);

    return GBM_OK;
}


int size_of_vector
(
    VEC_VEC_CATEGORIES &vec,
    int i
)
{
    return vec[i].size();
}



