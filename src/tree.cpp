//  GBM by Greg Ridgeway  Copyright (C) 2003
#include "tree.h"

CCARTTree::CCARTTree()
{
    pRootNode = NULL;
    pNodeFactory = NULL;
    dShrink = 1.0;
}


CCARTTree::~CCARTTree()
{
    if(pRootNode != NULL)
    {
        pRootNode->RecycleSelf(pNodeFactory);
    }
}


GBMRESULT CCARTTree::Initialize
(
    CNodeFactory *pNodeFactory
)
{
    GBMRESULT hr = GBM_OK;

    this->pNodeFactory = pNodeFactory;

    return hr;
}

    
GBMRESULT CCARTTree::Reset()
{
    GBMRESULT hr = GBM_OK;

    if(pRootNode != NULL)
    {
        // delete the old tree and start over
        hr = pRootNode->RecycleSelf(pNodeFactory);
    }
    if(GBM_FAILED(hr))
    {
        goto Error;
    }

    iBestNode = 0;
    dBestNodeImprovement = 0.0;

    schWhichNode = 0;

    pNewSplitNode    = NULL;
    pNewLeftNode     = NULL;
    pNewRightNode    = NULL;
    pNewMissingNode  = NULL;
    pInitialRootNode = NULL;

Cleanup:
    return hr;
Error:
    goto Cleanup;
}



//------------------------------------------------------------------------------
// Grows a regression tree
//------------------------------------------------------------------------------
GBMRESULT CCARTTree::grow
(
    double *adZ, 
    CDataset *pData, 
    double *adW,
    double *adF,
    unsigned long nTrain,
    unsigned long nBagged,
    double dLambda,
    unsigned long cMaxDepth,
    unsigned long cMinObsInNode,
    bool *afInBag,
    unsigned long *aiNodeAssign,
    CNodeSearch *aNodeSearch,
    VEC_P_NODETERMINAL &vecpTermNodes
)
{
    GBMRESULT hr = GBM_OK;

    #ifdef NOISY_DEBUG
    Rprintf("Growing tree\n");
    #endif

    if((adZ==NULL) || (pData==NULL) || (adW==NULL) || (adF==NULL) || 
       (cMaxDepth < 1))
    {
        hr = GBM_INVALIDARG;
        goto Error;
    }

    dSumZ = 0.0;
    dSumZ2 = 0.0;
    dTotalW = 0.0;

    #ifdef NOISY_DEBUG
    Rprintf("initial tree calcs\n");
    #endif
    for(iObs=0; iObs<nTrain; iObs++)
    {
        // aiNodeAssign tracks to which node each training obs belongs
        aiNodeAssign[iObs] = 0;
        if(afInBag[iObs])
        {
            // get the initial sums and sum of squares and total weight
            dSumZ += adW[iObs]*adZ[iObs];
            dSumZ2 += adW[iObs]*adZ[iObs]*adZ[iObs];
            dTotalW += adW[iObs];
        }
    }
    dError = dSumZ2-dSumZ*dSumZ/dTotalW;

    pInitialRootNode = pNodeFactory->GetNewNodeTerminal();
    pInitialRootNode->dPrediction = dSumZ/dTotalW;
    pInitialRootNode->dTrainW = dTotalW;
    vecpTermNodes.resize(2*cMaxDepth + 1,NULL); // accounts for missing nodes
    vecpTermNodes[0] = pInitialRootNode;
    pRootNode = pInitialRootNode;

    aNodeSearch[0].Set(dSumZ,dTotalW,nBagged,
                       pInitialRootNode,
                       &pRootNode,
                       pNodeFactory);

    // build the tree structure
    #ifdef NOISY_DEBUG
    Rprintf("Building tree 1 ");
    #endif
    cTotalNodeCount = 1;
    cTerminalNodes = 1;
    for(cDepth=0; cDepth<cMaxDepth; cDepth++)
    {
        #ifdef NOISY_DEBUG
        Rprintf("%d ",cDepth);
        #endif
        hr = GetBestSplit(pData,
                          nTrain,
                          aNodeSearch,
                          cTerminalNodes,
                          aiNodeAssign,
                          afInBag,
                          adZ,
                          adW,
                          iBestNode,
                          dBestNodeImprovement);
        if(GBM_FAILED(hr))
        {
            goto Error;
        }

        if(dBestNodeImprovement == 0.0)
        {
            break;
        }

        // setup the new nodes and add them to the tree
        hr = aNodeSearch[iBestNode].SetupNewNodes(pNewSplitNode,
                                                  pNewLeftNode,
                                                  pNewRightNode,
                                                  pNewMissingNode);
        cTotalNodeCount += 3;
        cTerminalNodes += 2;
        vecpTermNodes[iBestNode] = pNewLeftNode;
        vecpTermNodes[cTerminalNodes-2] = pNewRightNode;
        vecpTermNodes[cTerminalNodes-1] = pNewMissingNode;

        // assign observations to the correct node
        for(iObs=0; iObs < nTrain; iObs++)
        {
            iWhichNode = aiNodeAssign[iObs];
            if(iWhichNode==iBestNode)
            {
                schWhichNode = pNewSplitNode->WhichNode(pData,iObs);
                if(schWhichNode == 1) // goes right
                {
                    aiNodeAssign[iObs] = cTerminalNodes-2;
                }
                else if(schWhichNode == 0) // is missing
                {
                    aiNodeAssign[iObs] = cTerminalNodes-1;
                }
                // those to the left stay with the same node assignment
            }
        }

        // set up the node search for the new right node
        aNodeSearch[cTerminalNodes-2].Set(aNodeSearch[iBestNode].dBestRightSumZ,
                                          aNodeSearch[iBestNode].dBestRightTotalW,
                                          aNodeSearch[iBestNode].cBestRightN,
                                          pNewRightNode,
                                          &(pNewSplitNode->pRightNode),
                                          pNodeFactory);
        // set up the node search for the new missing node
        aNodeSearch[cTerminalNodes-1].Set(aNodeSearch[iBestNode].dBestMissingSumZ,
                                          aNodeSearch[iBestNode].dBestMissingTotalW,
                                          aNodeSearch[iBestNode].cBestMissingN,
                                          pNewMissingNode,
                                          &(pNewSplitNode->pMissingNode),
                                          pNodeFactory);
        // set up the node search for the new left node
        // must be done second since we need info for right node first
        aNodeSearch[iBestNode].Set(aNodeSearch[iBestNode].dBestLeftSumZ,
                                   aNodeSearch[iBestNode].dBestLeftTotalW,
                                   aNodeSearch[iBestNode].cBestLeftN,
                                   pNewLeftNode,
                                   &(pNewSplitNode->pLeftNode),
                                   pNodeFactory);

    } // end tree growing

    // DEBUG
    // Print();

Cleanup:
    return hr;
Error:
    goto Cleanup;
}


GBMRESULT CCARTTree::GetBestSplit
(
    CDataset *pData,
    unsigned long nTrain,
    CNodeSearch *aNodeSearch,
    unsigned long cTerminalNodes,
    unsigned long *aiNodeAssign,
    bool *afInBag,
    double *adZ,
    double *adW,
    unsigned long &iBestNode,
    double &dBestNodeImprovement
)
{
    GBMRESULT hr = GBM_OK;

    int iVar = 0;
    unsigned long iNode = 0;
    unsigned long iOrderObs = 0;
    unsigned long iWhichObs = 0;
    unsigned long cVarClasses = 0;
    double dX = 0.0;

    for(iVar=0; iVar < pData->cCols; iVar++)
    {
        cVarClasses = pData->acVarClasses[iVar];

        for(iNode=0; iNode < cTerminalNodes; iNode++)
        {
            hr = aNodeSearch[iNode].ResetForNewVar(iVar,cVarClasses);
        }

        // distribute the observations in order to the correct node search
        for(iOrderObs=0; iOrderObs < nTrain; iOrderObs++)
        {
            iWhichObs = pData->aiXOrder[iVar*nTrain + iOrderObs];
            if(afInBag[iWhichObs])
            {
                iNode = aiNodeAssign[iWhichObs];
                dX = pData->adX[iVar*(pData->cRows) + iWhichObs];
                hr = aNodeSearch[iNode].IncorporateObs
                     (dX,
                      adZ[iWhichObs],
                      adW[iWhichObs],
                      pData->alMonotoneVar[iVar]);
                if(GBM_FAILED(hr))
                {
                    goto Error;
                }
            }
        }
        for(iNode=0; iNode<cTerminalNodes; iNode++)
        {
            if(cVarClasses != 0) // evaluate if categorical split
            {
                hr = aNodeSearch[iNode].EvaluateCategoricalSplit();
            }
            aNodeSearch[iNode].WrapUpCurrentVariable();
        }
    }

    // search for the best split
    iBestNode = 0;
    dBestNodeImprovement = 0.0;
    for(iNode=0; iNode<cTerminalNodes; iNode++)
    {
        aNodeSearch[iNode].SetToSplit();
        if(aNodeSearch[iNode].BestImprovement() > dBestNodeImprovement)
        {
            iBestNode = iNode;
            dBestNodeImprovement = aNodeSearch[iNode].BestImprovement();
        }
    }

Cleanup:
    return hr;
Error:
    goto Cleanup;
}


GBMRESULT CCARTTree::GetNodeCount
(
    int &cNodes
)
{
    cNodes = cTotalNodeCount;
    return GBM_OK;
}



GBMRESULT CCARTTree::PredictValid
(
    CDataset *pData, 
    unsigned long nValid, 
    double *adFadj
)
{
    GBMRESULT hr = GBM_OK;
    int i=0;

    for(i=pData->cRows - nValid; i<pData->cRows; i++)
    {
        pRootNode->Predict(pData, i, adFadj[i]);
        adFadj[i] *= dShrink;
    }

    return hr;
}



GBMRESULT CCARTTree::Predict
(
    double *adX,
    unsigned long cRow, 
    unsigned long cCol, 
    unsigned long iRow, 
    double &dFadj
)
{

    if(pRootNode != NULL)
    {
        pRootNode->Predict(adX,cRow,cCol,iRow,dFadj);
        dFadj *= dShrink;
    }
    else
    {
        dFadj = 0.0;
    }

    return GBM_OK;
}



GBMRESULT CCARTTree::Adjust
(
    unsigned long *aiNodeAssign,
    double *adFadj,
    unsigned long cTrain,
    VEC_P_NODETERMINAL &vecpTermNodes,
    unsigned long cMinObsInNode
)
{
    unsigned long hr = GBM_OK;
    unsigned long iObs = 0;
    
    hr = pRootNode->Adjust(cMinObsInNode);
    if(GBM_FAILED(hr))
    {
        goto Error;
    }

    // predict for the training observations
    for(iObs=0; iObs<cTrain; iObs++)
    {
        adFadj[iObs] = vecpTermNodes[aiNodeAssign[iObs]]->dPrediction;
    }

Cleanup:
    return hr;
Error:
    goto Cleanup;
}


GBMRESULT CCARTTree::Print()
{
    GBMRESULT hr = GBM_OK;

    if(pRootNode != NULL)
    {
        pRootNode->PrintSubtree(0);
        Rprintf("shrinkage: %f\n",dShrink);
        Rprintf("initial error: %f\n\n",dError);
    }

    return hr;
}



GBMRESULT CCARTTree::GetVarRelativeInfluence
(
    double *adRelInf
)
{
    GBMRESULT hr = GBM_OK;

    if(pRootNode != NULL)
    {
        hr = pRootNode->GetVarRelativeInfluence(adRelInf);
        if(GBM_FAILED(hr))
        {
            goto Error;
        }
    }

Cleanup:
    return hr;
Error:
    goto Cleanup;

}



GBMRESULT CCARTTree::TransferTreeToRList
(
    CDataset *pData,
    int *aiSplitVar,
    double *adSplitPoint,
    int *aiLeftNode,
    int *aiRightNode,
    int *aiMissingNode,
    double *adErrorReduction,
    double *adWeight,
    double *adPred,
    VEC_VEC_CATEGORIES &vecSplitCodes,
    int cCatSplitsOld,
    double dShrinkage
)
{
    GBMRESULT hr = GBM_OK;

    int iNodeID = 0;

    if(pRootNode != NULL)
    {
        hr = pRootNode->TransferTreeToRList(iNodeID,
                                            pData,
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
                                            dShrinkage);
    }
    else
    {
        hr = GBM_FAIL;
    }
    return hr;
}


