//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//
//  File:       node_search.cpp
//
//------------------------------------------------------------------------------

#include "node_search.h"

CNodeSearch::CNodeSearch()
    :k_cMaxClasses(1024)
{
    iBestSplitVar = 0;
    dBestSplitValue = 0.0;
    fIsSplit = false;

    dBestMissingTotalW = 0.0;
    dCurrentMissingTotalW = 0.0;
    dBestMissingSumZ = 0.0;
    dCurrentMissingSumZ = 0.0;

    adGroupSumZ = NULL;
    adGroupW = NULL;
    acGroupN = NULL;
    adGroupMean = NULL;
    aiCurrentCategory = NULL;
    aiBestCategory = NULL;

    iRank = UINT_MAX;
}


CNodeSearch::~CNodeSearch()
{
    if(adGroupSumZ != NULL)
    {
        delete [] adGroupSumZ;
        adGroupSumZ = NULL;
    }
    if(adGroupW != NULL)
    {
        delete [] adGroupW;
        adGroupW = NULL;
    }
    if(acGroupN != NULL)
    {
        delete [] acGroupN;
        acGroupN = NULL;
    }
    if(adGroupMean != NULL)
    {
        delete [] adGroupMean;
        adGroupMean = NULL;
    }
    if(aiCurrentCategory != NULL)
    {
        delete [] aiCurrentCategory;
        aiCurrentCategory = NULL;
    }
    if(aiBestCategory != NULL)
    {
        delete [] aiBestCategory;
        aiBestCategory = NULL;
    }
}


GBMRESULT CNodeSearch::Initialize
(
    unsigned long cMinObsInNode
)
{
    GBMRESULT hr = GBM_OK;

    adGroupSumZ = new double[k_cMaxClasses];
    if(adGroupSumZ == NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }
    adGroupW = new double[k_cMaxClasses];
    if(adGroupW == NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }
    acGroupN = new ULONG[k_cMaxClasses];
    if(acGroupN == NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }
    adGroupMean = new double[k_cMaxClasses];
    if(adGroupMean == NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }
    aiCurrentCategory = new int[k_cMaxClasses];
    if(aiCurrentCategory == NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }
    aiBestCategory = new ULONG[k_cMaxClasses];
    if(aiBestCategory == NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }

    this->cMinObsInNode = cMinObsInNode;

Cleanup:
    return hr;
Error:
    goto Cleanup;
}


GBMRESULT CNodeSearch::IncorporateObs
(
    double dX,
    double dZ,
    double dW,
    long lMonotone
)
{
    GBMRESULT hr = GBM_OK;
    static double dWZ = 0.0;

    if(fIsSplit) goto Cleanup;

    dWZ = dW*dZ;

    if(ISNA(dX))
    {
        dCurrentMissingSumZ += dWZ;
        dCurrentMissingTotalW += dW;
        cCurrentMissingN++;
        dCurrentRightSumZ -= dWZ;
        dCurrentRightTotalW -= dW;
        cCurrentRightN--;
    }
    else if(cCurrentVarClasses == 0)   // variable is continuous
    {
        if(dLastXValue > dX)
        {
            error("Observations are not in order. gbm() was unable to build an index for the design matrix. Could be a bug in gbm or an unusual data type in data.\n");
            hr = GBM_FAIL;
            goto Error;
        }

        // Evaluate the current split
        // the newest observation is still in the right child
        dCurrentSplitValue = 0.5*(dLastXValue + dX);
        if((dLastXValue != dX) && 
            (cCurrentLeftN >= cMinObsInNode) && 
            (cCurrentRightN >= cMinObsInNode) &&
            ((lMonotone==0) ||
            (lMonotone*(dCurrentRightSumZ*dCurrentLeftTotalW - 
                        dCurrentLeftSumZ*dCurrentRightTotalW) > 0)))
        {
            dCurrentImprovement = 
                CNode::Improvement(dCurrentLeftTotalW,dCurrentRightTotalW,
                                    dCurrentMissingTotalW,
                                    dCurrentLeftSumZ,dCurrentRightSumZ,
                                    dCurrentMissingSumZ);
            if(dCurrentImprovement > dBestImprovement)
            {
                iBestSplitVar = iCurrentSplitVar;
                dBestSplitValue = dCurrentSplitValue;
                cBestVarClasses = 0;

                dBestLeftSumZ    = dCurrentLeftSumZ;
                dBestLeftTotalW  = dCurrentLeftTotalW;
                cBestLeftN       = cCurrentLeftN;
                dBestRightSumZ   = dCurrentRightSumZ;
                dBestRightTotalW = dCurrentRightTotalW;
                cBestRightN      = cCurrentRightN;
                dBestImprovement = dCurrentImprovement;
            }
        }

        // now move the new observation to the left
        // if another observation arrives we will evaluate this
        dCurrentLeftSumZ += dWZ;
        dCurrentLeftTotalW += dW;
        cCurrentLeftN++;
        dCurrentRightSumZ -= dWZ;
        dCurrentRightTotalW -= dW;
        cCurrentRightN--;

        dLastXValue = dX;
    }
    else // variable is categorical, evaluates later
    {
        adGroupSumZ[(unsigned long)dX] += dWZ;
        adGroupW[(unsigned long)dX] += dW;
        acGroupN[(unsigned long)dX] ++;
    }

Cleanup:
    return hr;
Error:
    goto Cleanup;
}



GBMRESULT CNodeSearch::Set
(
    double dSumZ,
    double dTotalW,
    unsigned long cTotalN,
    CNodeTerminal *pThisNode,
    CNode **ppParentPointerToThisNode,
    CNodeFactory *pNodeFactory
)
{
    GBMRESULT hr = GBM_OK;

    dInitSumZ = dSumZ;
    dInitTotalW = dTotalW;
    cInitN = cTotalN;

    dBestLeftSumZ       = 0.0;
    dBestLeftTotalW     = 0.0;
    cBestLeftN          = 0;
    dCurrentLeftSumZ    = 0.0;
    dCurrentLeftTotalW  = 0.0;
    cCurrentLeftN       = 0;

    dBestRightSumZ      = dSumZ;
    dBestRightTotalW    = dTotalW;
    cBestRightN         = cTotalN;
    dCurrentRightSumZ   = 0.0;
    dCurrentRightTotalW = dTotalW;
    cCurrentRightN      = cTotalN;

    dBestMissingSumZ      = 0.0;
    dBestMissingTotalW    = 0.0;
    cBestMissingN         = 0;
    dCurrentMissingSumZ   = 0.0;
    dCurrentMissingTotalW = 0.0;
    cCurrentMissingN      = 0;

    dBestImprovement    = 0.0;
    iBestSplitVar       = UINT_MAX;

    dCurrentImprovement = 0.0;
    iCurrentSplitVar    = UINT_MAX;
    dCurrentSplitValue  = -HUGE_VAL;

    fIsSplit = false;

    this->pThisNode = pThisNode;
    this->ppParentPointerToThisNode = ppParentPointerToThisNode;
    this->pNodeFactory = pNodeFactory;

    return hr;
}


GBMRESULT CNodeSearch::ResetForNewVar
(
    unsigned long iWhichVar,
    long cCurrentVarClasses
)
{
    GBMRESULT hr = GBM_OK;
    long i=0;

    if(fIsSplit) goto Cleanup;

    for(i=0; i<cCurrentVarClasses; i++)
    {
        adGroupSumZ[i] = 0.0;
        adGroupW[i] = 0.0;
        acGroupN[i] = 0;
    }

    iCurrentSplitVar = iWhichVar;
    this->cCurrentVarClasses = cCurrentVarClasses;

    dCurrentLeftSumZ      = 0.0;
    dCurrentLeftTotalW    = 0.0;
    cCurrentLeftN         = 0;
    dCurrentRightSumZ     = dInitSumZ;
    dCurrentRightTotalW   = dInitTotalW;
    cCurrentRightN        = cInitN;
    dCurrentMissingSumZ   = 0.0;
    dCurrentMissingTotalW = 0.0;
    cCurrentMissingN      = 0;

    dCurrentImprovement = 0.0;

    dLastXValue = -HUGE_VAL;

Cleanup:
    return hr;
}



GBMRESULT CNodeSearch::WrapUpCurrentVariable()
{
    GBMRESULT hr = GBM_OK;
    if(iCurrentSplitVar == iBestSplitVar)
    {
        if(cCurrentMissingN > 0)
        {
            dBestMissingSumZ   = dCurrentMissingSumZ;
            dBestMissingTotalW = dCurrentMissingTotalW;
            cBestMissingN      = cCurrentMissingN;
        }
        else // DEBUG: consider a weighted average with parent node?
        {
            dBestMissingSumZ   = dInitSumZ;
            dBestMissingTotalW = dInitTotalW;
            cBestMissingN      = 0;
        }
    }

    return hr;
}



GBMRESULT CNodeSearch::EvaluateCategoricalSplit()
{
    GBMRESULT hr = GBM_OK;
    long i=0;
    long j=0;
    unsigned long cFiniteMeans = 0;

    if(fIsSplit) goto Cleanup;

    if(cCurrentVarClasses == 0)
    {
        hr = GBM_INVALIDARG;
        goto Error;
    }

    cFiniteMeans = 0;
    for(i=0; i<cCurrentVarClasses; i++)
    {
        aiCurrentCategory[i] = i;
        if(adGroupW[i] != 0.0)
        {
            adGroupMean[i] = adGroupSumZ[i]/adGroupW[i];
            cFiniteMeans++;
        }
        else
        {
            adGroupMean[i] = HUGE_VAL;
        }
    }

    rsort_with_index(adGroupMean,aiCurrentCategory,cCurrentVarClasses);

    // if only one group has a finite mean it will not consider
    // might be all are missing so no categories enter here
    for(i=0; (cFiniteMeans>1) && ((ULONG)i<cFiniteMeans-1); i++)
    {
        dCurrentSplitValue = (double)i;

        dCurrentLeftSumZ    += adGroupSumZ[aiCurrentCategory[i]];
        dCurrentLeftTotalW  += adGroupW[aiCurrentCategory[i]];
        cCurrentLeftN       += acGroupN[aiCurrentCategory[i]];
        dCurrentRightSumZ   -= adGroupSumZ[aiCurrentCategory[i]];
        dCurrentRightTotalW -= adGroupW[aiCurrentCategory[i]];
        cCurrentRightN      -= acGroupN[aiCurrentCategory[i]];

        dCurrentImprovement = 
            CNode::Improvement(dCurrentLeftTotalW,dCurrentRightTotalW,
                               dCurrentMissingTotalW,
                               dCurrentLeftSumZ,dCurrentRightSumZ,
                               dCurrentMissingSumZ);
        if((cCurrentLeftN >= cMinObsInNode) && 
           (cCurrentRightN >= cMinObsInNode) &&
           (dCurrentImprovement > dBestImprovement))
        {
            dBestSplitValue = dCurrentSplitValue;
            if(iBestSplitVar != iCurrentSplitVar)
            {
                iBestSplitVar = iCurrentSplitVar;
                cBestVarClasses = cCurrentVarClasses;
                for(j=0; j<cCurrentVarClasses; j++)
                {
                    aiBestCategory[j] = aiCurrentCategory[j];
                }
            }

            dBestLeftSumZ      = dCurrentLeftSumZ;
            dBestLeftTotalW    = dCurrentLeftTotalW;
            cBestLeftN         = cCurrentLeftN;
            dBestRightSumZ     = dCurrentRightSumZ;
            dBestRightTotalW   = dCurrentRightTotalW;
            cBestRightN        = cCurrentRightN;
            dBestImprovement   = dCurrentImprovement;
        }
    }

Cleanup:
    return hr;
Error:
    goto Cleanup;
}




GBMRESULT CNodeSearch::SetupNewNodes
(
    PCNodeNonterminal &pNewSplitNode,
    PCNodeTerminal &pNewLeftNode,
    PCNodeTerminal &pNewRightNode,
    PCNodeTerminal &pNewMissingNode
)
{
    GBMRESULT hr = GBM_OK;
    CNodeContinuous *pNewNodeContinuous = NULL;
    CNodeCategorical *pNewNodeCategorical = NULL;
    unsigned long i=0;

    pNewLeftNode    = pNodeFactory->GetNewNodeTerminal();
    pNewRightNode   = pNodeFactory->GetNewNodeTerminal();
    pNewMissingNode = pNodeFactory->GetNewNodeTerminal();

    // set up a continuous split
    if(cBestVarClasses==0)
    {
        pNewNodeContinuous = pNodeFactory->GetNewNodeContinuous();

        pNewNodeContinuous->dSplitValue = dBestSplitValue;
        pNewNodeContinuous->iSplitVar = iBestSplitVar;

        pNewSplitNode = pNewNodeContinuous;
    }
    else
    {
        // get a new categorical node and its branches
        pNewNodeCategorical = pNodeFactory->GetNewNodeCategorical();

        // set up the categorical split
        pNewNodeCategorical->iSplitVar = iBestSplitVar;
        pNewNodeCategorical->cLeftCategory = (ULONG)dBestSplitValue + 1;
        pNewNodeCategorical->aiLeftCategory = 
            new ULONG[pNewNodeCategorical->cLeftCategory];
        for(i=0; i<pNewNodeCategorical->cLeftCategory; i++)
        {
            pNewNodeCategorical->aiLeftCategory[i] = aiBestCategory[i];
        }

        pNewSplitNode = pNewNodeCategorical;
    }

    *ppParentPointerToThisNode = pNewSplitNode;

    pNewSplitNode->dPrediction  = pThisNode->dPrediction;
    pNewSplitNode->dImprovement = dBestImprovement;
    pNewSplitNode->dTrainW      = pThisNode->dTrainW;
    pNewSplitNode->pLeftNode    = pNewLeftNode;
    pNewSplitNode->pRightNode   = pNewRightNode;
    pNewSplitNode->pMissingNode = pNewMissingNode;

    pNewLeftNode->dPrediction    = dBestLeftSumZ/dBestLeftTotalW;
    pNewLeftNode->dTrainW        = dBestLeftTotalW;
    pNewLeftNode->cN             = cBestLeftN;
    pNewRightNode->dPrediction   = dBestRightSumZ/dBestRightTotalW;
    pNewRightNode->dTrainW       = dBestRightTotalW;
    pNewRightNode->cN            = cBestRightN;
    pNewMissingNode->dPrediction = dBestMissingSumZ/dBestMissingTotalW;
    pNewMissingNode->dTrainW     = dBestMissingTotalW;
    pNewMissingNode->cN          = cBestMissingN;

    pThisNode->RecycleSelf(pNodeFactory);

    return hr;
}
