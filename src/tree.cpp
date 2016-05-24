//  GBM by Greg Ridgeway  Copyright (C) 2003
//-----------------------------------
// Includes
//-----------------------------------
#include <algorithm>
#include "tree.h"

//----------------------------------------
// Function Members - Public
//----------------------------------------
CCARTTree::CCARTTree(double shrinkage, long depth):
kTreeDepth_(depth), kShrinkage_(shrinkage)
{
    rootnode_ = NULL;
    totalnodecount_ = 1;
    error_ = 0.0;

    // Calculate original
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
}



//------------------------------------------------------------------------------
// Grows a regression tree
//------------------------------------------------------------------------------
void CCARTTree::grow
(
 double* residuals,
 const CDataset& kData,
 const double* kFuncEstimates,
 unsigned long min_num_node_obs,
 std::vector<unsigned long>& data_node_assigns,
 CNodeSearch& nodesearcher
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
		data_node_assigns[obs_num] = 0;

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
  nodesearcher.set_search_rootnode(*rootnode_);

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
    nodesearcher.GenerateAllSplits(terminalnode_ptrs_, kData, &(residuals[0]), data_node_assigns);
    double bestImprov = nodesearcher.CalcImprovementAndSplit(terminalnode_ptrs_, kData, data_node_assigns);

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
 const std::vector<unsigned long>& kDataNodeAssigns,
 double* delta_estimates,
 unsigned long min_num_node_obs
)
{
	unsigned long obs_num = 0;

	rootnode_->Adjust(min_num_node_obs);

	// predict for the training observations
	for(obs_num=0; obs_num<kDataNodeAssigns.size(); obs_num++)
	{
		delta_estimates[obs_num] = terminalnode_ptrs_[kDataNodeAssigns[obs_num]]->prediction;
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


CNode* CCARTTree::GetRootNode()
{
	return rootnode_;
}

const CNode* CCARTTree::GetRootNode() const
{
	return rootnode_;
}


