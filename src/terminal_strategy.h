//------------------------------------------------------------------------------
//
//  File:       terminalStrategy.h
//
//  Description: Strategies for terminal nodes
//
//	Author: 	James Hickey
//------------------------------------------------------------------------------

#ifndef TERMINALSTRATEGY_H
#define TERMINALSTRATEGY_H

//------------------------------
// Includes
//------------------------------
#include "dataset.h"
#include "generic_node_strategy.h"
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
	TerminalStrategy(CNode* node):nodecontext_(node){};

	//---------------------
	// Public destructor
	//---------------------
	~TerminalStrategy(){nodecontext_ = NULL;};

	//---------------------
	// Public Functions
	//---------------------
	void Adjust(unsigned long cMinObsInNode){return;};
	void Predict(const CDataset &data,
		    unsigned long iRow,
		    double &dFadj)
	{
		dFadj = nodecontext_->prediction;
	}
	void GetVarRelativeInfluence(double* adRelInf){return;}
	void PrintSubTree(unsigned long Indent)
	{
		  for(unsigned long i=0; i< Indent; i++) Rprintf("  ");
		  Rprintf("N=%f, Prediction=%f *\n",
			  nodecontext_->totalweight,
			  nodecontext_->prediction);
	}
	signed char WhichNode(const CDataset& data, unsigned long iObs)
	{
		signed char ReturnValue = 0;
		double dX = data.x_value(iObs, nodecontext_->split_var);
		 if(!ISNA(dX))
			{
				if(dX < nodecontext_->splitvalue)
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
		adSplitPoint[iNodeID] = dShrinkage*nodecontext_->prediction;
		aiLeftNode[iNodeID] = -1;
		aiRightNode[iNodeID] = -1;
		aiMissingNode[iNodeID] = -1;
		adErrorReduction[iNodeID] = 0.0;
		adWeight[iNodeID] = nodecontext_->totalweight;
		adPred[iNodeID] = dShrinkage*nodecontext_->prediction;

		iNodeID++;
	}

private:
	CNode* nodecontext_;
};
#endif // TERMINALSTRATEGY_H
