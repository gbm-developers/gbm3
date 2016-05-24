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
 double *adZ,
 const CDataset& data,
 const double *adF,
 unsigned long cMinObsInNode,
 std::vector<unsigned long>& aiNodeAssign,
 CNodeSearch& aNodeSearch
)
{
#ifdef NOISY_DEBUG
  Rprintf("Growing tree\n");
#endif

	if((adZ==NULL) || (data.weight_ptr()==NULL) || (adF==NULL) ||
	 (kTreeDepth_ < 1))
	{
	  throw GBM::InvalidArgument();
	}

  double dSumZ = 0.0;
  double dSumZ2 = 0.0;
  double dTotalW = 0.0;

#ifdef NOISY_DEBUG
  Rprintf("initial tree calcs\n");
#endif

  	  // Move to data -- FOR TIME BEING
	for(unsigned long iObs=0; iObs<data.get_trainsize(); iObs++)
	{
		// aiNodeAssign tracks to which node each training obs belongs
		aiNodeAssign[iObs] = 0;

		if(data.get_bag_element(iObs))
		{
			// get the initial sums and sum of squares and total weight
			dSumZ += data.weight_ptr()[iObs]*adZ[iObs];
			dSumZ2 += data.weight_ptr()[iObs]*adZ[iObs]*adZ[iObs];
			dTotalW += data.weight_ptr()[iObs];
		}
	}

  error_ = dSumZ2-dSumZ*dSumZ/dTotalW;
  rootnode_ = new CNode(NodeDef(dSumZ, dTotalW, data.get_total_in_bag()));
  terminalnode_ptrs_[0] = rootnode_;
  aNodeSearch.set_search_rootnode(*rootnode_);

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
    aNodeSearch.GenerateAllSplits(terminalnode_ptrs_, data, &(adZ[0]), aiNodeAssign);
    double bestImprov = aNodeSearch.CalcImprovementAndSplit(terminalnode_ptrs_, data, aiNodeAssign);

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
 const CDataset &data,
 unsigned long nValid,
 double *adFadj
 )
{
  unsigned int i=0;
  for(i=data.nrow() - nValid; i<data.nrow(); i++)
    {
      rootnode_->Predict(data, i, adFadj[i]);
      adFadj[i] *= kShrinkage_;
    }
}

void CCARTTree::Adjust
(
 const std::vector<unsigned long>& aiNodeAssign,
 double *adFadj,
 unsigned long cMinObsInNode
)
{
	unsigned long iObs = 0;

	rootnode_->Adjust(cMinObsInNode);

	// predict for the training observations
	for(iObs=0; iObs<aiNodeAssign.size(); iObs++)
	{
		adFadj[iObs] = terminalnode_ptrs_[aiNodeAssign[iObs]]->prediction;
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


