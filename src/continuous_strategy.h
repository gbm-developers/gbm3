//------------------------------------------------------------------------------
//
//  File:       continuousStrategy.h
//
//  Description: strategies for continuous splits.
//
//	Author: 	James Hickey
//------------------------------------------------------------------------------

#ifndef CONTINUOUSSTRATEGY_H
#define CONTINUOUSSTRATEGY_H

//------------------------------
// Includes
//------------------------------
#include "dataset.h"
#include "node.h"
#include "generic_node_strategy.h"
#include <Rcpp.h>


//------------------------------
// Class Definition
//------------------------------
class ContinuousStrategy:public GenericNodeStrategy
{
public:
	//----------------------
	// Public Constructors
	//----------------------
	ContinuousStrategy(CNode* node):node_context_(node){};

	//---------------------
	// Public destructor
	//---------------------
	~ContinuousStrategy(){node_context_=NULL;};

	//---------------------
	// Public Functions
	//---------------------
	void Adjust(unsigned long cMinObsInNode)
	{
		node_context_->left_node_ptr->Adjust(cMinObsInNode);
		node_context_->right_node_ptr->Adjust(cMinObsInNode);

		if((node_context_->missing_node_ptr->splittype == none) && (node_context_->missing_node_ptr->numobs < cMinObsInNode))
		{
			node_context_->prediction = ((node_context_->left_node_ptr->totalweight)*(node_context_->left_node_ptr->prediction) +
				 (node_context_->right_node_ptr->totalweight)*(node_context_->right_node_ptr->prediction))/
			(node_context_->left_node_ptr->totalweight + node_context_->right_node_ptr->totalweight);
			node_context_->missing_node_ptr->prediction = node_context_->prediction;
		}
		else
		{
			node_context_->missing_node_ptr->Adjust(cMinObsInNode);
			node_context_->prediction =
			((node_context_->left_node_ptr->totalweight)*(node_context_->left_node_ptr->prediction) +
			(node_context_->right_node_ptr->totalweight)*  (node_context_->right_node_ptr->prediction) +
			(node_context_->missing_node_ptr->totalweight)*(node_context_->missing_node_ptr->prediction))/
			(node_context_->left_node_ptr->totalweight + node_context_->right_node_ptr->totalweight + node_context_->missing_node_ptr->totalweight);
		}
	}
	void Predict(const CDataset &data,
			    unsigned long iRow,
			    double &dFadj)
	{
		signed char schWhichNode = ContinuousStrategy::WhichNode(data,iRow);
		if(schWhichNode == -1)
		{
		  node_context_->left_node_ptr->Predict(data, iRow, dFadj);
		}
		else if(schWhichNode == 1)
		{
		  node_context_->right_node_ptr->Predict(data, iRow, dFadj);
		}
		else
		{
		  node_context_->missing_node_ptr->Predict(data, iRow, dFadj);
		}
	}
	void GetVarRelativeInfluence(double* adRelInf)
	{
		adRelInf[node_context_->split_var] += node_context_->improvement;
		node_context_->left_node_ptr->GetVarRelativeInfluence(adRelInf);
		node_context_->right_node_ptr->GetVarRelativeInfluence(adRelInf);
	}
	void PrintSubTree(unsigned long Indent)
	{
		const std::size_t cLeftCategory = node_context_->leftcategory.size();

		for(unsigned long i=0; i< Indent; i++) Rprintf("  ");
		Rprintf("N=%f, Improvement=%f, Prediction=%f, NA pred=%f\n",
		  node_context_->totalweight,
		  node_context_->improvement,
		  node_context_->prediction,
		  (node_context_->missing_node_ptr == NULL ? 0.0 : node_context_->missing_node_ptr->prediction));

		for(unsigned long i=0; i< Indent; i++) Rprintf("  ");
		Rprintf("V%d in ",node_context_->split_var);
		for(unsigned long i=0; i<cLeftCategory; i++)
		  {
			Rprintf("%d", node_context_->leftcategory[i]);
			if(i<cLeftCategory-1) Rprintf(",");
		  }
		Rprintf("\n");
		node_context_->left_node_ptr->PrintSubtree((Indent+1));

		for(unsigned long i=0; i< Indent; i++) Rprintf("  ");
		Rprintf("V%d not in ", node_context_->split_var);
		for(unsigned long i=0; i<cLeftCategory; i++)
		  {
			Rprintf("%d", node_context_->leftcategory[i]);
			if(i<cLeftCategory-1) Rprintf(",");
		  }
		Rprintf("\n");
		node_context_->right_node_ptr->PrintSubtree(Indent+1);

		for(unsigned long i=0; i< Indent; i++) Rprintf("  ");
		Rprintf("missing\n");
		node_context_->missing_node_ptr->PrintSubtree(Indent+1);
	}
	signed char WhichNode(const CDataset& data, unsigned long iObs)
	{
		signed char ReturnValue = 0;
		double dX = data.x_value(iObs, node_context_->split_var);
		 if(!ISNA(dX))
			{
				if(dX < node_context_->splitvalue)
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
		aiSplitVar[iThisNodeID] = node_context_->split_var;
		adSplitPoint[iThisNodeID] = node_context_->splitvalue;
		adErrorReduction[iThisNodeID] = node_context_->improvement;
		adWeight[iThisNodeID] = node_context_->totalweight;
		adPred[iThisNodeID] = dShrinkage*node_context_->prediction;

		iNodeID++;
		aiLeftNode[iThisNodeID] = iNodeID;
		node_context_->left_node_ptr->TransferTreeToRList(iNodeID,
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
		node_context_->right_node_ptr->TransferTreeToRList(iNodeID,
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
		node_context_->missing_node_ptr->TransferTreeToRList(iNodeID,
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
	CNode* node_context_;
};
#endif // CONTINUOUSSTRATEGY_H

