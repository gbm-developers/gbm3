//------------------------------------------------------------------------------
//
//  File:       terminalStrategy.h
//
//  Description: Strategies for terminal nodes
//
//	Author: 	James Hickey
//------------------------------------------------------------------------------

#ifndef __terminalStrategy_h__
#define __terminalStrategy_h__

//------------------------------
// Includes
//------------------------------
#include "dataset.h"
#include "genericNodeStrategy.h"
#include "node.h"
#include <Rcpp.h>


//------------------------------
// Class Definition
//------------------------------
class TerminalStrategy: public GenericNodeStrategy
{
public:
	//----------------------
	// Public Constructors
	//----------------------
	TerminalStrategy(CNode* node):nodeContext(node){};

	//---------------------
	// Public destructor
	//---------------------
	~TerminalStrategy(){nodeContext = NULL;};

	//---------------------
	// Public Functions
	//---------------------
	void Adjust(unsigned long cMinObsInNode){return;};
	void Predict(const CDataset &data,
		    unsigned long iRow,
		    double &dFadj)
	{
		dFadj = nodeContext->dPrediction;
	}
	void GetVarRelativeInfluence(double* adRelInf){return;}
	void PrintSubTree(unsigned long Indent)
	{
		  for(long i=0; i< Indent; i++) Rprintf("  ");
		  Rprintf("N=%f, Prediction=%f *\n",
			  nodeContext->dTrainW,
			  nodeContext->dPrediction);
	}
	signed char WhichNode(const CDataset& data, unsigned long iObs)
	{
		signed char ReturnValue = 0;
		double dX = data.x_value(iObs, nodeContext->iSplitVar);
		 if(!ISNA(dX))
			{
				if(dX < nodeContext->dSplitValue)
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
		aiSplitVar[iNodeID] = -1;
		adSplitPoint[iNodeID] = dShrinkage*nodeContext->dPrediction;
		aiLeftNode[iNodeID] = -1;
		aiRightNode[iNodeID] = -1;
		aiMissingNode[iNodeID] = -1;
		adErrorReduction[iNodeID] = 0.0;
		adWeight[iNodeID] = nodeContext->dTrainW;
		adPred[iNodeID] = dShrinkage*nodeContext->dPrediction;

		iNodeID++;
	}

private:
	CNode* nodeContext;
};
#endif //__terminalStrategy_h__
