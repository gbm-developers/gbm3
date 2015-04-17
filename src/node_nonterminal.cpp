//  GBM by Greg Ridgeway  Copyright (C) 2003

#include "node_nonterminal.h"

CNodeNonterminal::CNodeNonterminal()
{
    pLeftNode = NULL;
    pRightNode = NULL;
    iSplitVar = 0;
    dImprovement = 0.0;
    pMissingNode = NULL;
}


CNodeNonterminal::~CNodeNonterminal()
{

}


void CNodeNonterminal::Adjust
(
    unsigned long cMinObsInNode
)
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



void CNodeNonterminal::Predict
(
    CDataset &data,
    unsigned long iRow,
    double &dFadj
)
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


void CNodeNonterminal::Predict
(
    double *adX,
    unsigned long cRow,
    unsigned long cCol,
    unsigned long iRow,
    double &dFadj
)
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


void CNodeNonterminal::GetVarRelativeInfluence
(
    double *adRelInf
)
{
  adRelInf[iSplitVar] += dImprovement;
  pLeftNode->GetVarRelativeInfluence(adRelInf);
  pRightNode->GetVarRelativeInfluence(adRelInf);
}


