//  GBM by Greg Ridgeway  Copyright (C) 2003

#include "node.h"

CNode::CNode()
{
    dPrediction = 0.0;
    dTrainW = 0.0;
    isTerminal = false;
    cN = 0;

	pLeftNode = NULL;
	pRightNode = NULL;
	iSplitVar = 0;
	dImprovement = 0.0;
	pMissingNode = NULL;
}


CNode::~CNode()
{
	// Each node is responsible for deleting its
	// children
    delete pLeftNode;
    delete pRightNode;
    delete pMissingNode;
}


void CNode::Adjust
(
    unsigned long cMinObsInNode
)
{
	// Only adjust if node is not terminal
	if(!isTerminal)
	{
		pLeftNode->Adjust(cMinObsInNode);
		pRightNode->Adjust(cMinObsInNode);

		if(pMissingNode->isTerminal && (pMissingNode->cN < cMinObsInNode))
		{
			dPrediction = ((pLeftNode->dTrainW)*(pLeftNode->dPrediction) +
				 (pRightNode->dTrainW)*(pRightNode->dPrediction))/
			(pLeftNode->dTrainW + pRightNode->dTrainW);
			pMissingNode->dPrediction = dPrediction;
		}
		else
		{
			pMissingNode->Adjust(cMinObsInNode);
			dPrediction =
			((pLeftNode->dTrainW)*   (pLeftNode->dPrediction) +
			(pRightNode->dTrainW)*  (pRightNode->dPrediction) +
			(pMissingNode->dTrainW)*(pMissingNode->dPrediction))/
			(pLeftNode->dTrainW + pRightNode->dTrainW + pMissingNode->dTrainW);
		}
	}

}



void CNode::Predict
(
    const CDataset &data,
    unsigned long iRow,
    double &dFadj
)
{
	// If node is terminal set the function adjustment to the
	// prediction.  Else move down the tree.
	if(isTerminal)
	{
		dFadj = dPrediction;
	}
	else
	{
		signed char schWhichNode = WhichNode(data,iRow);
		if(schWhichNode == -1)
		{
		  pLeftNode->Predict(data, iRow, dFadj);
		}
		else if(schWhichNode == 1)
		{
		  pRightNode->Predict(data, iRow, dFadj);
		}
		else
		{
		  pMissingNode->Predict(data, iRow, dFadj);
		}
	}

}


void CNode::Predict
(
    double *adX,
    unsigned long cRow,
    unsigned long cCol,
    unsigned long iRow,
    double &dFadj
)
{
	// If node is terminal set the function adjustment to the
	// prediction.  Else move down the tree.
	if(isTerminal)
	{
		dFadj = dPrediction;
	}
	else
	{
		signed char schWhichNode = WhichNode(adX,cRow,cCol,iRow);
		if(schWhichNode == -1)
		{
		  pLeftNode->Predict(adX,cRow,cCol,iRow,dFadj);
		}
		else if(schWhichNode == 1)
		{
		  pRightNode->Predict(adX,cRow,cCol,iRow,dFadj);
		}
		else
		{
		  pMissingNode->Predict(adX,cRow,cCol,iRow,dFadj);
		}
	}

}


void CNode::GetVarRelativeInfluence
(
    double *adRelInf
)
{
	// Relative influence of split variable only updated in non-terminal nodes
	if(!isTerminal)
	{
		adRelInf[iSplitVar] += dImprovement;
		pLeftNode->GetVarRelativeInfluence(adRelInf);
		pRightNode->GetVarRelativeInfluence(adRelInf);
	}

}

void CNode::ApplyShrinkage(double dLambda)
{
	dPrediction *= dLambda;
}




