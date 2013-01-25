//  GBM by Greg Ridgeway  Copyright (C) 2003

#include "node.h"

CNode::CNode()
{
    dPrediction = 0.0;
    dTrainW = 0.0;
    isTerminal = false;
}


CNode::~CNode()
{
    // the nodes get deleted by deleting the node factory
}


GBMRESULT CNode::Adjust
(
    unsigned long cMinObsInNode
)
{
    GBMRESULT hr = GBM_NOTIMPL;
    return hr;
}


GBMRESULT CNode::Predict
(
    CDataset *pData, 
    unsigned long iRow, 
    double &dFadj
)
{
    GBMRESULT hr = GBM_NOTIMPL;
    return hr;
}


double CNode::TotalError()
{
    GBMRESULT hr = GBM_NOTIMPL;
    return hr;
}


GBMRESULT CNode::PrintSubtree
(
    unsigned long cIndent
)
{
    GBMRESULT hr = GBM_NOTIMPL;
    return hr;
}


GBMRESULT CNode::GetVarRelativeInfluence
(
    double *adRelInf
)
{
    GBMRESULT hr = GBM_NOTIMPL;
    return hr;
}


GBMRESULT CNode::TransferTreeToRList
(
    int &iNodeID,
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
    return GBM_NOTIMPL;
}


