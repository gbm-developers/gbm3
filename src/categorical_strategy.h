//------------------------------------------------------------------------------
//
//  File:       categoricalStrategy.h
//
//  Description: strategy for categorical splits.
//
//	Author: 	James Hickey
//------------------------------------------------------------------------------

#ifndef CATEGORICALSTRATEGY_H
#define CATEGORICALSTRATEGY_H

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
class CategoricalStrategy:public GenericNodeStrategy
{
public:
	//----------------------
	// Public Constructors
	//----------------------
	CategoricalStrategy(CNode* node):node_context_(node){};

	//---------------------
	// Public destructor
	//---------------------
	~CategoricalStrategy(){node_context_=NULL;};

	//---------------------
	// Public Functions
	//---------------------
	void Adjust(unsigned long min_num_node_obs)
	{
		node_context_->left_node_ptr->Adjust(min_num_node_obs);
		node_context_->right_node_ptr->Adjust(min_num_node_obs);

		if((node_context_->missing_node_ptr->splittype == kNone) && (node_context_->missing_node_ptr->numobs < min_num_node_obs))
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
				    unsigned long rownum,
				    double &deltaestimate)
	{
		signed char whichnode = CategoricalStrategy::WhichNode(kData,rownum);
		if(whichnode == -1)
		{
		  node_context_->left_node_ptr->Predict(kData, rownum, deltaestimate);
		}
		else if(whichnode == 1)
		{
		  node_context_->right_node_ptr->Predict(kData, rownum, deltaestimate);
		}
		else
		{
		  node_context_->missing_node_ptr->Predict(kData, rownum, deltaestimate);
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
		signed char returnvalue = 0;
		double xval = kData.x_value(obs_num, node_context_->split_var);

    	if(!ISNA(xval))
		{
		  if(std::find(node_context_->leftcategory.begin(),
			   node_context_->leftcategory.end(),
			   (unsigned long)xval) != node_context_->leftcategory.end())
			{
				returnvalue = -1;
			}
			else
			{
				returnvalue = 1;
			}
		}
		// if missing value returns 0

		return returnvalue;
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
		VecOfVectorCategories &splitcodes_vec,
		int prev_categorical_splits,
		double shrinkage
	)
	{
		int thisnode_id = nodeid;
		unsigned long cat_splits = splitcodes_vec.size();
		unsigned long i = 0;
		int levels = kData.varclass(node_context_->split_var);
		const std::size_t leftcategory = node_context_->leftcategory.size();

		splitvar[thisnode_id] = node_context_->split_var;
		splitvalues[thisnode_id] = cat_splits+prev_categorical_splits; // 0 based
		error_reduction[thisnode_id] = node_context_->improvement;
		weights[thisnode_id] = node_context_->totalweight;
		predictions[thisnode_id] = shrinkage*node_context_->prediction;

		splitcodes_vec.push_back(VectorCategories());

		splitcodes_vec[cat_splits].resize(levels,1);
		for(i=0; i<leftcategory; i++)
		  {
			splitcodes_vec[cat_splits][node_context_->leftcategory[i]] = -1;
		  }

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
#endif // CATEGORICALSTRATEGY_H
