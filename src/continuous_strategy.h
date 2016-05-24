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
	void Adjust(unsigned long min_num_node_obs)
	{
		node_context_->left_node_ptr->Adjust(min_num_node_obs);
		node_context_->right_node_ptr->Adjust(min_num_node_obs);

		if((node_context_->missing_node_ptr->splittype == none) && (node_context_->missing_node_ptr->numobs < min_num_node_obs))
		{
			node_context_->prediction = ((node_context_->left_node_ptr->totalweight)*(node_context_->left_node_ptr->prediction) +
				 (node_context_->right_node_ptr->totalweight)*(node_context_->right_node_ptr->prediction))/
			(node_context_->left_node_ptr->totalweight + node_context_->right_node_ptr->totalweight);
			node_context_->missing_node_ptr->prediction = node_context_->prediction;
		}
		else
		{
			node_context_->missing_node_ptr->Adjust(min_num_node_obs);
			node_context_->prediction =
			((node_context_->left_node_ptr->totalweight)*(node_context_->left_node_ptr->prediction) +
			(node_context_->right_node_ptr->totalweight)*  (node_context_->right_node_ptr->prediction) +
			(node_context_->missing_node_ptr->totalweight)*(node_context_->missing_node_ptr->prediction))/
			(node_context_->left_node_ptr->totalweight + node_context_->right_node_ptr->totalweight + node_context_->missing_node_ptr->totalweight);
		}
	}
	void Predict(const CDataset &kData,
			    unsigned long row_num,
			    double &deltafunc_est)
	{
		signed char whichnode = ContinuousStrategy::WhichNode(kData,row_num);
		if(whichnode == -1)
		{
		  node_context_->left_node_ptr->Predict(kData, row_num, deltafunc_est);
		}
		else if(whichnode == 1)
		{
		  node_context_->right_node_ptr->Predict(kData, row_num, deltafunc_est);
		}
		else
		{
		  node_context_->missing_node_ptr->Predict(kData, row_num, deltafunc_est);
		}
	}
	void GetVarRelativeInfluence(double* relative_influence)
	{
		relative_influence[node_context_->split_var] += node_context_->improvement;
		node_context_->left_node_ptr->GetVarRelativeInfluence(relative_influence);
		node_context_->right_node_ptr->GetVarRelativeInfluence(relative_influence);
	}
	void PrintSubTree(unsigned long indent)
	{
		const std::size_t leftcategory = node_context_->leftcategory.size();

		for(unsigned long i=0; i< indent; i++) Rprintf("  ");
		Rprintf("N=%f, Improvement=%f, Prediction=%f, NA pred=%f\n",
		  node_context_->totalweight,
		  node_context_->improvement,
		  node_context_->prediction,
		  (node_context_->missing_node_ptr == NULL ? 0.0 : node_context_->missing_node_ptr->prediction));

		for(unsigned long i=0; i< indent; i++) Rprintf("  ");
		Rprintf("V%d in ",node_context_->split_var);
		for(unsigned long i=0; i<leftcategory; i++)
		  {
			Rprintf("%d", node_context_->leftcategory[i]);
			if(i<leftcategory-1) Rprintf(",");
		  }
		Rprintf("\n");
		node_context_->left_node_ptr->PrintSubtree((indent+1));

		for(unsigned long i=0; i< indent; i++) Rprintf("  ");
		Rprintf("V%d not in ", node_context_->split_var);
		for(unsigned long i=0; i<leftcategory; i++)
		  {
			Rprintf("%d", node_context_->leftcategory[i]);
			if(i<leftcategory-1) Rprintf(",");
		  }
		Rprintf("\n");
		node_context_->right_node_ptr->PrintSubtree(indent+1);

		for(unsigned long i=0; i< indent; i++) Rprintf("  ");
		Rprintf("missing\n");
		node_context_->missing_node_ptr->PrintSubtree(indent+1);
	}
	signed char WhichNode(const CDataset& kData, unsigned long obs_num)
	{
		signed char return_value = 0;
		double xval = kData.x_value(obs_num, node_context_->split_var);
		 if(!ISNA(xval))
			{
				if(xval < node_context_->splitvalue)
				{
					return_value = -1;
				}
				else
				{
					return_value = 1;
				}
			}
			// if missing value returns 0

			return return_value;
	}
	void TransferTreeToRList
	(
		int &nodeid,
		const CDataset &kData,
		int* splitvar,
		double* splitvalues,
		int* leftnodes,
		int* rightnodes,
		int* missingnodes,
		double* error_reduction,
		double* weights,
		double* predictions,
		VEC_VEC_CATEGORIES &splitcodes_vec,
		int prev_categorical_splits,
		double shrinkage
	)
	{
		int thisnode_id = nodeid;
		splitvar[thisnode_id] = node_context_->split_var;
		splitvalues[thisnode_id] = node_context_->splitvalue;
		error_reduction[thisnode_id] = node_context_->improvement;
		weights[thisnode_id] = node_context_->totalweight;
		predictions[thisnode_id] = shrinkage*node_context_->prediction;

		nodeid++;
		leftnodes[thisnode_id] = nodeid;
		node_context_->left_node_ptr->TransferTreeToRList(nodeid,
					 kData,
					 splitvar,
					 splitvalues,
					 leftnodes,
					 rightnodes,
					 missingnodes,
					 error_reduction,
					 weights,
					 predictions,
					 splitcodes_vec,
					 prev_categorical_splits,
					 shrinkage);

		rightnodes[thisnode_id] = nodeid;
		node_context_->right_node_ptr->TransferTreeToRList(nodeid,
					  kData,
					  splitvar,
					  splitvalues,
					  leftnodes,
					  rightnodes,
					  missingnodes,
					  error_reduction,
					  weights,
					  predictions,
					  splitcodes_vec,
					  prev_categorical_splits,
					  shrinkage);
		missingnodes[thisnode_id] = nodeid;
		node_context_->missing_node_ptr->TransferTreeToRList(nodeid,
						kData,
						splitvar,
						splitvalues,
						leftnodes,
						rightnodes,
						missingnodes,
						error_reduction,
						weights,
						predictions,
						splitcodes_vec,
						prev_categorical_splits,
						shrinkage);
	}

private:
	CNode* node_context_;
};
#endif // CONTINUOUSSTRATEGY_H

