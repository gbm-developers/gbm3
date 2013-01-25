//  GBM by Greg Ridgeway  Copyright (C) 2003

#include "node_factory.h"

CNodeFactory::CNodeFactory()
{
}


CNodeFactory::~CNodeFactory()
{
    #ifdef NOISY_DEBUG
    Rprintf("destructing node factory\n");
    #endif
}


GBMRESULT CNodeFactory::Initialize
(
    unsigned long cDepth
)
{
    GBMRESULT hr = GBM_OK;
    unsigned long i = 0;

    for(i=0; i<NODEFACTORY_NODGBM_RESERVE; i++)
    {
        TerminalStack.push(&(aBlockTerminal[i]));
        ContinuousStack.push(&(aBlockContinuous[i]));
        CategoricalStack.push(&(aBlockCategorical[i]));
    }

    return hr;
}


CNodeTerminal* CNodeFactory::GetNewNodeTerminal()
{
    if(TerminalStack.empty())
    {
        #ifdef NOISY_DEBUG
        Rprintf("Terminal stack is empty\n");
        #endif
        pNodeTerminalTemp = NULL;
    }
    else
    {
        pNodeTerminalTemp = TerminalStack.top();
        TerminalStack.pop();

        pNodeTerminalTemp->dPrediction = 0.0;
    }
    return pNodeTerminalTemp;
}


CNodeContinuous* CNodeFactory::GetNewNodeContinuous()
{
    if(ContinuousStack.empty())
    {
        #ifdef NOISY_DEBUG
        Rprintf("Continuous stack is empty\n");
        #endif
        pNodeContinuousTemp = NULL;
    }
    else
    {
        pNodeContinuousTemp = ContinuousStack.top(); 
        ContinuousStack.pop();

        pNodeContinuousTemp->dPrediction = 0.0;
        pNodeContinuousTemp->dImprovement = 0.0;
        pNodeContinuousTemp->pMissingNode = NULL;
        pNodeContinuousTemp->pLeftNode = NULL;
        pNodeContinuousTemp->pRightNode = NULL;
        pNodeContinuousTemp->iSplitVar = 0;
        pNodeContinuousTemp->dSplitValue = 0.0;
    }

    return pNodeContinuousTemp;
}


CNodeCategorical* CNodeFactory::GetNewNodeCategorical()
{
    if(CategoricalStack.empty())
    {
        #ifdef NOISY_DEBUG
        Rprintf("Categorical stack is empty\n");
        #endif
        pNodeCategoricalTemp = NULL;
    }
    else
    {
        pNodeCategoricalTemp = CategoricalStack.top();
        CategoricalStack.pop();

        pNodeCategoricalTemp->dPrediction = 0.0;
        pNodeCategoricalTemp->dImprovement = 0.0;
        pNodeCategoricalTemp->pMissingNode = NULL;
        pNodeCategoricalTemp->pLeftNode = NULL;
        pNodeCategoricalTemp->pRightNode = NULL;
        pNodeCategoricalTemp->iSplitVar = 0;
        pNodeCategoricalTemp->aiLeftCategory = NULL;
        pNodeCategoricalTemp->cLeftCategory = 0;
    }

    return pNodeCategoricalTemp;
}


GBMRESULT CNodeFactory::RecycleNode
(
    CNodeTerminal *pNode
)
{
    if(pNode != NULL)
    {
        TerminalStack.push(pNode);
    }
    return GBM_OK;
}

GBMRESULT CNodeFactory::RecycleNode
(
    CNodeContinuous *pNode
)
{
    if(pNode != NULL)
    {
        if(pNode->pLeftNode != NULL) pNode->pLeftNode->RecycleSelf(this);
        if(pNode->pRightNode != NULL) pNode->pRightNode->RecycleSelf(this);
        if(pNode->pMissingNode != NULL) pNode->pMissingNode->RecycleSelf(this);
        ContinuousStack.push(pNode);
    }
    return GBM_OK;
}

GBMRESULT CNodeFactory::RecycleNode
(
    CNodeCategorical *pNode
)
{
    if(pNode != NULL)
    {
        if(pNode->pLeftNode != NULL) pNode->pLeftNode->RecycleSelf(this);
        if(pNode->pRightNode != NULL) pNode->pRightNode->RecycleSelf(this);
        if(pNode->pMissingNode != NULL) pNode->pMissingNode->RecycleSelf(this);
        if(pNode->aiLeftCategory != NULL)
        {
            delete [] pNode->aiLeftCategory;
            pNode->aiLeftCategory = NULL;
        }
        CategoricalStack.push(pNode);
    }

    return GBM_OK;
}


