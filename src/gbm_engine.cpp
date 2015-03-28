//  GBM by Greg Ridgeway  Copyright (C) 2003
//#define NOISY_DEBUG
#include <algorithm>

#include "gbm_engine.h"

CGBM::CGBM()
{
    afInBag = NULL;
    aNodeSearch = NULL;

    cDepth = 0;
    cMinObsInNode = 0;
    dBagFraction = 0.0;
    dLambda = 0.0;
    fInitialized = false;
    cTotalInBag = 0;
    cTrain = 0;
    cFeatures = 0;
    cValid = 0;

    pData = NULL;
    pDist = NULL;
    pNodeFactory = NULL;
    ptreeTemp = NULL;
}


CGBM::~CGBM()
{
    delete[] afInBag;
    delete[] aNodeSearch;
    delete ptreeTemp;
    delete pNodeFactory;
}


GBMRESULT CGBM::Initialize
(
    CDataset *pData,
    CDistribution *pDist,
    double dLambda,
    unsigned long cTrain,
    unsigned long cFeatures,
    double dBagFraction,
    unsigned long cDepth,
    unsigned long cMinObsInNode,
    unsigned long cNumClasses,
    int cGroups
)
{
    GBMRESULT hr = GBM_OK;
    unsigned long i=0;

    if(pData == NULL)
    {
        hr = GBM_INVALIDARG;
        goto Error;
    }
    if(pDist == NULL)
    {
        hr = GBM_INVALIDARG;
        goto Error;
    }

    this->pData = pData;
    this->pDist = pDist;
    this->dLambda = dLambda;
    this->cTrain = cTrain;
    this->cFeatures = cFeatures;
    this->dBagFraction = dBagFraction;
    this->cDepth = cDepth;
    this->cMinObsInNode = cMinObsInNode;
    this->cGroups = cGroups;

    // allocate the tree structure
    ptreeTemp = new CCARTTree;

    cValid = pData->cRows - cTrain;
    cTotalInBag = (unsigned long)(dBagFraction*cTrain);
    adZ.assign((pData->cRows) * cNumClasses, 0);
    adFadj.assign((pData->cRows) * cNumClasses, 0);

    pNodeFactory = new CNodeFactory();
    pNodeFactory->Initialize(cDepth);
    ptreeTemp->Initialize(pNodeFactory);

    // array for flagging those observations in the bag
    afInBag = new bool[cTrain];

    // aiNodeAssign tracks to which node each training obs belongs
    aiNodeAssign.resize(cTrain);
    // NodeSearch objects help decide which nodes to split
    aNodeSearch = new CNodeSearch[2*cDepth+1];

    for(i=0; i<2*cDepth+1; i++)
    {
        aNodeSearch[i].Initialize(cMinObsInNode);
    }
    vecpTermNodes.resize(2*cDepth+1,NULL);

    fInitialized = true;

Cleanup:
    return hr;
Error:
    goto Cleanup;
}




GBMRESULT CGBM::Predict
(
    unsigned long iVar,
    unsigned long cTrees,
    double *adF,
    double *adX,
    unsigned long cLength
)
{
    GBMRESULT hr = GBM_OK;


    return hr;
}


GBMRESULT CGBM::Predict
(
    double *adX,
    unsigned long cRow,
    unsigned long cCol,
    unsigned long cTrees,
    double *adF
)
{
    GBMRESULT hr = GBM_OK;


    return hr;
}



GBMRESULT CGBM::GetVarRelativeInfluence
(
    double *adRelInf,
    unsigned long cTrees
)
{
    GBMRESULT hr = GBM_OK;

    std::fill(adRelInf, adRelInf + pData->cCols, 0);

    return hr;
}


GBMRESULT CGBM::PrintTree()
{
    GBMRESULT hr = GBM_OK;

    hr = ptreeTemp->Print();
    if(GBM_FAILED(hr)) goto Error;

Cleanup:
    return hr;
Error:
    goto Cleanup;
}


GBMRESULT CGBM::iterate
(
    double *adF,
    double &dTrainError,
    double &dValidError,
    double &dOOBagImprove,
    int &cNodes,
    int cNumClasses,
    int cClassIdx
)
{
    GBMRESULT hr = GBM_OK;
    unsigned long i = 0;
    unsigned long cBagged = 0;
    int cIdxOff = cClassIdx * (cTrain + cValid);

 //   for(i=0; i < cTrain + cIdxOff; i++){ adF[i] = 0;}
    if(!fInitialized)
    {
        hr = GBM_FAIL;
        goto Error;
    }

    dTrainError = 0.0;
    dValidError = 0.0;
    dOOBagImprove = 0.0;

    vecpTermNodes.assign(2*cDepth+1,NULL);

    // randomly assign observations to the Bag

    if (cClassIdx == 0)
    {
        if (!IsPairwise())
        {
            // regular instance based training
            for(i=0; i<cTrain; i++) /* && (cBagged < cTotalInBag); i++) */
            {
                if(unif_rand()*(cTrain-i) < cTotalInBag-cBagged)
                {
                    afInBag[i] = true;
                    cBagged++;
                }
                else
                {
                    afInBag[i] = false;
                }
/*                if (cBagged >= cTotalInBag){
                    break;
                } */
            }
	    std::fill(afInBag + i, afInBag + cTrain, false);
        }
        else
        {
            // for pairwise training, sampling is per group
            // therefore, we will not have exactly cTotalInBag instances

            double dLastGroup = -1;
            bool chosen = false;
            unsigned int cBaggedGroups = 0;
            unsigned int cSeenGroups   = 0;
            unsigned int cTotalGroupsInBag = (unsigned long)(dBagFraction * cGroups);
            if (cTotalGroupsInBag <= 0)
            {
                cTotalGroupsInBag = 1;
            }
            for(i=0; i<cTrain; i++)
            {
                const double dGroup = pData->adMisc[i];
                if (dGroup != dLastGroup)
                {
                    if (cBaggedGroups >= cTotalGroupsInBag)
                    {
                        break;
                    }

                    // Group changed, make a new decision
                    chosen = (unif_rand()*(cGroups - cSeenGroups) < cTotalGroupsInBag - cBaggedGroups);
                    if (chosen)
                    {
                        cBaggedGroups++;
                    }
                    dLastGroup = dGroup;
                    cSeenGroups++;
                }
                if (chosen)
                {
                    afInBag[i] = true;
                    cBagged++;
                }
                else
                {
                    afInBag[i] = false;
                }
            }
            // the remainder is not in the bag
	    std::fill(afInBag + i, afInBag + cTrain, false);
        }
    }

#ifdef NOISY_DEBUG
    Rprintf("Compute working response\n");
#endif

    hr = pDist->ComputeWorkingResponse(pData->adY,
                                       pData->adMisc,
                                       pData->adOffset,
                                       adF,
                                       &adZ[0],
                                       pData->adWeight,
                                       afInBag,
                                       cTrain,
                                       cIdxOff);

    if(GBM_FAILED(hr))
    {
        goto Error;
    }

#ifdef NOISY_DEBUG
    Rprintf("Reset tree\n");
#endif
    ptreeTemp->Reset();
#ifdef NOISY_DEBUG
    Rprintf("grow tree\n");
#endif

    hr = ptreeTemp->grow(&(adZ[cIdxOff]), pData, &(pData->adWeight[cIdxOff]),
                         &(adFadj[cIdxOff]), cTrain, cFeatures, cTotalInBag, dLambda, cDepth,
                         cMinObsInNode, afInBag, aiNodeAssign, aNodeSearch,
                         vecpTermNodes);

    if(GBM_FAILED(hr))
    {
        goto Error;
    }

#ifdef NOISY_DEBUG
    Rprintf("get node count\n");
#endif
    hr = ptreeTemp->GetNodeCount(cNodes);
    if(GBM_FAILED(hr))
    {
        goto Error;
    }

    // Now I have adF, adZ, and vecpTermNodes (new node assignments)
    // Fit the best constant within each terminal node
#ifdef NOISY_DEBUG
    Rprintf("fit best constant\n");
#endif

    hr = pDist->FitBestConstant(pData->adY,
                                pData->adMisc,
                                pData->adOffset,
                                pData->adWeight,
                                &adF[0],
                                &adZ[0],
                                aiNodeAssign,
                                cTrain,
                                vecpTermNodes,
                                (2*cNodes+1)/3, // number of terminal nodes
                                cMinObsInNode,
                                afInBag,
                                &adFadj[0],
                                cIdxOff);

    if(GBM_FAILED(hr))
    {
        goto Error;
    }

    // update training predictions
    // fill in missing nodes where N < cMinObsInNode
    hr = ptreeTemp->Adjust(aiNodeAssign,&(adFadj[cIdxOff]),cTrain,
                           vecpTermNodes,cMinObsInNode);
    if(GBM_FAILED(hr))
    {
        goto Error;
    }
    ptreeTemp->SetShrinkage(dLambda);

    if (cClassIdx == (cNumClasses - 1))
    {
        dOOBagImprove = pDist->BagImprovement(pData->adY,
                                              pData->adMisc,
                                              pData->adOffset,
                                              pData->adWeight,
                                              &adF[0],
                                              &adFadj[0],
                                              afInBag,
                                              dLambda,
                                              cTrain);
    }

    // update the training predictions
    for(i=0; i < cTrain; i++)
    {
        int iIdx = i + cIdxOff;
        adF[iIdx] += dLambda * adFadj[iIdx];
    }

    dTrainError = pDist->Deviance(pData->adY,
                                  pData->adMisc,
                                  pData->adOffset,
                                  pData->adWeight,
                                  adF,
                                  cTrain,
                                  cIdxOff);

    // update the validation predictions
    hr = ptreeTemp->PredictValid(pData,cValid,&(adFadj[cIdxOff]));

    for(i=cTrain; i < cTrain+cValid; i++)
    {
        adF[i + cIdxOff] += adFadj[i + cIdxOff];
    }

    if(pData->fHasOffset)
    {
        dValidError =
            pDist->Deviance(pData->adY,
                            pData->adMisc,
                            pData->adOffset,
                            pData->adWeight,
                            adF,
                            cValid,
                            cIdxOff + cTrain);
    }
    else
    {
        dValidError = pDist->Deviance(pData->adY,
                                      pData->adMisc,
                                      NULL,
                                      pData->adWeight,
                                      adF,
                                      cValid,
                                      cIdxOff + cTrain);
    }

Cleanup:
    return hr;
Error:
    goto Cleanup;
}


GBMRESULT CGBM::TransferTreeToRList
(
    int *aiSplitVar,
    double *adSplitPoint,
    int *aiLeftNode,
    int *aiRightNode,
    int *aiMissingNode,
    double *adErrorReduction,
    double *adWeight,
    double *adPred,
    VEC_VEC_CATEGORIES &vecSplitCodes,
    int cCatSplitsOld
)
{
    GBMRESULT hr = GBM_OK;

    hr = ptreeTemp->TransferTreeToRList(pData,
                                        aiSplitVar,
                                        adSplitPoint,
                                        aiLeftNode,
                                        aiRightNode,
                                        aiMissingNode,
                                        adErrorReduction,
                                        adWeight,
                                        adPred,
                                        vecSplitCodes,
                                        cCatSplitsOld,
                                        dLambda);

    return hr;
}


