//------------------------------------------------------------------------------
//
//  File:       categoricalStrategy.h
//
//  Description: strategy for categorical splits.
//
//	Author: 	James Hickey
//------------------------------------------------------------------------------

#ifndef __categoricalStrategy_h__
#define __categoricalStrategy_h__

//------------------------------
// Includes
//------------------------------
#include "dataset.h"
#include "node.h"
#include "genericNodeStrategy.h"
#include <Rcpp.h>


//------------------------------
// Class Definition
//------------------------------
class CategoricalStrategy:public GenericNodeStrategy
{
public:
	//----------------------
	// Public Constructors
	//----------------------
	CategoricalStrategy(CNode* node):nodeContext(node){};

	//---------------------
	// Public destructor
	//---------------------
	~CategoricalStrategy(){nodeContext=NULL;};

	//---------------------
	// Public Functions
	//---------------------
	void Adjust(unsigned long cMinObsInNode)
	{
		nodeContext->pLeftNode->Adjust(cMinObsInNode);
		nodeContext->pRightNode->Adjust(cMinObsInNode);

		if((nodeContext->pMissingNode->splitType == none) && (nodeContext->pMissingNode->cN < cMinObsInNode))
		{
			nodeContext->dPrediction = ((nodeContext->pLeftNode->dTrainW)*(nodeContext->pLeftNode->dPrediction) +
				 (nodeContext->pRightNode->dTrainW)*(nodeContext->pRightNode->dPrediction))/
			(nodeContext->pLeftNode->dTrainW + nodeContext->pRightNode->dTrainW);
			nodeContext->pMissingNode->dPrediction = nodeContext->dPrediction;
		}
		else
		{
			nodeContext->pMissingNode->Adjust(cMinObsInNode);
			nodeContext->dPrediction =
			((nodeContext->pLeftNode->dTrainW)*(nodeContext->pLeftNode->dPrediction) +
			(nodeContext->pRightNode->dTrainW)*  (nodeContext->pRightNode->dPrediction) +
			(nodeContext->pMissingNode->dTrainW)*(nodeContext->pMissingNode->dPrediction))/
			(nodeContext->pLeftNode->dTrainW + nodeContext->pRightNode->dTrainW + nodeContext->pMissingNode->dTrainW);
		}
	}
	void Predict(const CDataset &data,
				    unsigned long iRow,
				    double &dFadj)
	{
		signed char schWhichNode = CategoricalStrategy::WhichNode(data,iRow);
		if(schWhichNode == -1)
		{
		  nodeContext->pLeftNode->Predict(data, iRow, dFadj);
		}
		else if(schWhichNode == 1)
		{
		  nodeContext->pRightNode->Predict(data, iRow, dFadj);
		}
		else
		{
		  nodeContext->pMissingNode->Predict(data, iRow, dFadj);
		}
	}
	void GetVarRelativeInfluence(double* adRelInf)
	{
		adRelInf[nodeContext->iSplitVar] += nodeContext->dImprovement;
		nodeContext->pLeftNode->GetVarRelativeInfluence(adRelInf);
		nodeContext->pRightNode->GetVarRelativeInfluence(adRelInf);
	}
	void PrintSubTree(unsigned long Indent)
	{
		const std::size_t cLeftCategory = nodeContext->aiLeftCategory.size();

		for(long i=0; i< Indent; i++) Rprintf("  ");
		Rprintf("N=%f, Improvement=%f, Prediction=%f, NA pred=%f\n",
		  nodeContext->dTrainW,
		  nodeContext->dImprovement,
		  nodeContext->dPrediction,
		  (nodeContext->pMissingNode == NULL ? 0.0 : nodeContext->pMissingNode->dPrediction));

		for(long i=0; i< Indent; i++) Rprintf("  ");
		Rprintf("V%d in ",nodeContext->iSplitVar);
		for(long i=0; i<cLeftCategory; i++)
		  {
			Rprintf("%d", nodeContext->aiLeftCategory[i]);
			if(i<cLeftCategory-1) Rprintf(",");
		  }
		Rprintf("\n");
		nodeContext->pLeftNode->PrintSubtree((Indent+1));

		for(long i=0; i< Indent; i++) Rprintf("  ");
		Rprintf("V%d not in ", nodeContext->iSplitVar);
		for(long i=0; i<cLeftCategory; i++)
		  {
			Rprintf("%d", nodeContext->aiLeftCategory[i]);
			if(i<cLeftCategory-1) Rprintf(",");
		  }
		Rprintf("\n");
		nodeContext->pRightNode->PrintSubtree(Indent+1);

		for(long i=0; i< Indent; i++) Rprintf("  ");
		Rprintf("missing\n");
		nodeContext->pMissingNode->PrintSubtree(Indent+1);
	}
	signed char WhichNode(const CDataset& data, unsigned long iObs)
	{
		signed char ReturnValue = 0;
		double dX = data.x_value(iObs, nodeContext->iSplitVar);

    	if(!ISNA(dX))
		{
		  if(std::find(nodeContext->aiLeftCategory.begin(),
			   nodeContext->aiLeftCategory.end(),
			   (ULONG)dX) != nodeContext->aiLeftCategory.end())
			{
				ReturnValue = -1;
			}
			else
			{
				ReturnValue = 1;
			}
		}
		// if missing value returns 0

		return ReturnValue;
	}
	void TransferTreeToRList
	(
		int &iNodeID,
		const CDataset &data,
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
		int iThisNodeID = iNodeID;
		unsigned long cCatSplits = vecSplitCodes.size();
		unsigned long i = 0;
		int cLevels = data.varclass(nodeContext->iSplitVar);
		const std::size_t cLeftCategory = nodeContext->aiLeftCategory.size();

		aiSplitVar[iThisNodeID] = nodeContext->iSplitVar;
		adSplitPoint[iThisNodeID] = cCatSplits+cCatSplitsOld; // 0 based
		adErrorReduction[iThisNodeID] = nodeContext->dImprovement;
		adWeight[iThisNodeID] = nodeContext->dTrainW;
		adPred[iThisNodeID] = dShrinkage*nodeContext->dPrediction;

		vecSplitCodes.push_back(VEC_CATEGORIES());

		vecSplitCodes[cCatSplits].resize(cLevels,1);
		for(i=0; i<cLeftCategory; i++)
		  {
			vecSplitCodes[cCatSplits][nodeContext->aiLeftCategory[i]] = -1;
		  }

		iNodeID++;
		aiLeftNode[iThisNodeID] = iNodeID;
		nodeContext->pLeftNode->TransferTreeToRList(iNodeID,
					 data,
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
		aiRightNode[iThisNodeID] = iNodeID;
		nodeContext->pRightNode->TransferTreeToRList(iNodeID,
					  data,
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

		aiMissingNode[iThisNodeID] = iNodeID;
		nodeContext->pMissingNode->TransferTreeToRList(iNodeID,
						data,
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

private:
	CNode* nodeContext;
};
#endif //__categoricalStrategy_h__
