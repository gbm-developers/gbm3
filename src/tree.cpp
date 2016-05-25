//  GBM by Greg Ridgeway  Copyright (C) 2003
//-----------------------------------
// Includes
//-----------------------------------
#include <algorithm>
#include "tree.h"

//----------------------------------------
// Function Members - Public
//----------------------------------------
CCARTTree::CCARTTree(TreeParams treeconfig) :
new_node_searcher_(treeconfig.depth, treeconfig.numberdatacolumns, treeconfig.min_obs_in_node),
kTreeDepth_(treeconfig.depth), kShrinkage_(treeconfig.shrinkage)
{
    rootnode_ = NULL;
    totalnodecount_ = 1;
    error_ = 0.0;
    min_num_node_obs_ = treeconfig.min_obs_in_node;
    data_node_assignment_.resize(treeconfig.num_trainrows, 0);
}


CCARTTree::~CCARTTree()
{
	delete rootnode_;
}

void CCARTTree::Reset()
{
  delete rootnode_;
  rootnode_ = NULL;
  terminalnode_ptrs_.resize(2*kTreeDepth_ + 1, NULL);
  totalnodecount_ = 1;
  new_node_searcher_.reset();

}



//------------------------------------------------------------------------------
// Grows a regression tree
//------------------------------------------------------------------------------
void CCARTTree::Grow
(
 double* residuals,
 const CDataset& kData,
 const double* kFuncEstimates
)
{
#ifdef NOISY_DEBUG
  Rprintf("Growing tree\n");
#endif

	if((residuals==NULL) || (kData.weight_ptr()==NULL) || (kFuncEstimates==NULL) ||
	 (kTreeDepth_ < 1))
	{
	  throw GBM::InvalidArgument();
	}

  double sumz = 0.0;
  double sum_zsquared = 0.0;
  double totalw = 0.0;

#ifdef NOISY_DEBUG
  Rprintf("initial tree calcs\n");
#endif

  	  // Move to data -- FOR TIME BEING
	for(unsigned long obs_num=0; obs_num<kData.get_trainsize(); obs_num++)
	{
		// aiNodeAssign tracks to which node each training obs belongs
		data_node_assignment_[obs_num] = 0;

		if(kData.get_bag_element(obs_num))
		{
			// get the initial sums and sum of squares and total weight
			sumz += kData.weight_ptr()[obs_num]*residuals[obs_num];
			sum_zsquared += kData.weight_ptr()[obs_num]*residuals[obs_num]*residuals[obs_num];
			totalw += kData.weight_ptr()[obs_num];
		}
	}

  error_ = sum_zsquared-sumz*sumz/totalw;
  rootnode_ = new CNode(NodeDef(sumz, totalw, kData.get_total_in_bag()));
  terminalnode_ptrs_[0] = rootnode_;
  new_node_searcher_.set_search_rootnode(*rootnode_);

  // build the tree structure
#ifdef NOISY_DEBUG
  Rprintf("Building tree 1 ");
#endif

  for(long cDepth=0; cDepth < kTreeDepth_; cDepth++)
  {
#ifdef NOISY_DEBUG
      Rprintf("%d ",cDepth);
#endif
      
    // Generate all splits
    new_node_searcher_.GenerateAllSplits(terminalnode_ptrs_, kData, &(residuals[0]), data_node_assignment_);
    double bestImprov = new_node_searcher_.CalcImprovementAndSplit(terminalnode_ptrs_, kData, data_node_assignment_);

    // Make the best split if possible
	if(bestImprov == 0.0)
	{
	  break;
	}
	// setup the new nodes and add them to the tree

	totalnodecount_ += 3;

  } // end tree growing

    // DEBUG
    // Print();
}

void CCARTTree::PredictValid
(
 const CDataset &kData,
 unsigned long num_validation_points,
 double* delta_estimates
 )
{
  unsigned int i=0;
  for(i=kData.nrow() - num_validation_points; i<kData.nrow(); i++)
    {
      rootnode_->Predict(kData, i, delta_estimates[i]);
      delta_estimates[i] *= kShrinkage_;
    }
}

void CCARTTree::Adjust
(
 double* delta_estimates
)
{
	rootnode_->Adjust(min_num_node_obs_);

	// predict for the training observations
	for(unsigned long obs_num=0; obs_num<data_node_assignment_.size(); obs_num++)
	{
		delta_estimates[obs_num] = terminalnode_ptrs_[data_node_assignment_[obs_num]]->prediction;
	}
}


void CCARTTree::Print()
{
    if(rootnode_)
    {
      rootnode_->PrintSubtree(0);
      Rprintf("shrinkage: %f\n",kShrinkage_);
      Rprintf("initial error: %f\n\n",error_);
    }
}

//-----------------------------------
// Function: TransferTreeToRList
//
// Returns: none
//
// Description:
//
// Parameters:
//
//-----------------------------------
void CCARTTree::TransferTreeToRList(const CDataset &kData,
	     int* splitvar,
	     double* splitvalues,
	     int* leftnodes,
	     int* rightnodes,
	     int* missingnodes,
	     double* error_reduction,
	     double* weights,
	     double* predictions,
	     VEC_VEC_CATEGORIES &splitcodes_vec,
	     int prev_categorical_splits)
{
	int nodeid = 0;
	if(rootnode_)
	{
		rootnode_->TransferTreeToRList(nodeid,
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
									   kShrinkage_);
	}
	else
	{
	  throw GBM::Failure("Can't transfer to list - RootNode does not exist.");
	}
}


