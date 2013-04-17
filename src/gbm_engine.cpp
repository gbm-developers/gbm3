//  GBM by Greg Ridgeway  Copyright (C) 2003
//#define NOISY_DEBUG
#include "gbm_engine.h"

CGBM::CGBM()
{
    adFadj = NULL;
    adZ = NULL;
    afInBag = NULL;
    aiNodeAssign = NULL;
    aNodeSearch = NULL;

    cDepth = 0;
    cMinObsInNode = 0;
    dBagFraction = 0.0;
    dLambda = 0.0;
    fInitialized = false;
    cTotalInBag = 0;
    cTrain = 0;
    cValid = 0;

    pData = NULL;
    pDist = NULL;
    pNodeFactory = NULL;
    ptreeTemp = NULL;
}


CGBM::~CGBM()
{
    if(adFadj != NULL)
    {
        delete [] adFadj;
        adFadj = NULL;
    }
    if(adZ != NULL)
    {
        delete [] adZ;
        adZ = NULL;
    }
    if(afInBag != NULL)
    {
        delete [] afInBag;
        afInBag = NULL;
    }
    if(aiNodeAssign != NULL)
    {
        delete [] aiNodeAssign;
        aiNodeAssign = NULL;
    }
    if(aNodeSearch != NULL)
    {
        delete [] aNodeSearch;
        aNodeSearch = NULL;
    }
    if(ptreeTemp != NULL)
    {
        delete ptreeTemp;
        ptreeTemp = NULL;
    }
    // must delete the node factory last!!! at least after deleting trees
    if(pNodeFactory != NULL)
    {
        delete pNodeFactory;
        pNodeFactory = NULL;
    }
}


GBMRESULT CGBM::Initialize
(
    CDataset *pData,
    CDistribution *pDist,
    double dLambda,
    unsigned long cTrain,
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
    this->dBagFraction = dBagFraction;
    this->cDepth = cDepth;
    this->cMinObsInNode = cMinObsInNode;
    this->cGroups = cGroups;

    // allocate the tree structure
    ptreeTemp = new CCARTTree;
    if(ptreeTemp == NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }

    cValid = pData->cRows - cTrain;
    cTotalInBag = (unsigned long)(dBagFraction*cTrain);
    adZ = new double[(pData->cRows) * cNumClasses];

    if(adZ == NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }
    adFadj = new double[(pData->cRows) * cNumClasses];
    if(adFadj == NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }

    for (i=0; i<(pData->cRows)*cNumClasses; i++)
    {
        adFadj[i] = 0.0;
    }

    pNodeFactory = new CNodeFactory();
    if(pNodeFactory == NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }
    hr = pNodeFactory->Initialize(cDepth);
    if(GBM_FAILED(hr))
    {
        goto Error;
    }
    ptreeTemp->Initialize(pNodeFactory);

    // array for flagging those observations in the bag
    afInBag = new bool[cTrain];
    if(afInBag==NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }
    // aiNodeAssign tracks to which node each training obs belongs
    aiNodeAssign = new ULONG[cTrain];
    if(aiNodeAssign==NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }
    // NodeSearch objects help decide which nodes to split
    aNodeSearch = new CNodeSearch[2*cDepth+1];
    if(aNodeSearch==NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }
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
    int iVar=0;

    for(iVar=0; iVar<pData->cCols; iVar++)
    {
        adRelInf[iVar] = 0.0;
    }

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
            // the remainder is not in the bag
            for( ; i<cTrain; i++)
            {
                afInBag[i] = false;
            }
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
            for( ; i<cTrain; i++)
            {
                afInBag[i] = false;
            }
        }
    }

#ifdef NOISY_DEBUG
    Rprintf("Compute working response\n");
#endif

    hr = pDist->ComputeWorkingResponse(pData->adY,
                                       pData->adMisc,
                                       pData->adOffset,
                                       adF,
                                       adZ,
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
    hr = ptreeTemp->Reset();
#ifdef NOISY_DEBUG
    Rprintf("grow tree\n");
#endif

    hr = ptreeTemp->grow(&(adZ[cIdxOff]), pData, &(pData->adWeight[cIdxOff]),
                         &(adFadj[cIdxOff]), cTrain, cTotalInBag, dLambda, cDepth,
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
                                adF,
                                adZ,
                                aiNodeAssign,
                                cTrain,
                                vecpTermNodes,
                                (2*cNodes+1)/3, // number of terminal nodes
                                cMinObsInNode,
                                afInBag,
                                adFadj,
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
                                              adF,
                                              adFadj,
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


