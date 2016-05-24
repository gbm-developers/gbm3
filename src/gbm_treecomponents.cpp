//-----------------------------------
//
// File: gbmTreeComps.cpp
//
// Description: class that implements the public methods of the
//    gbm engine tree components.
//
//-----------------------------------

//------------------------------
// Includes
//------------------------------
#include "gbm_treecomponents.h"

//----------------------------------------
// Function Members - Public
//----------------------------------------
//-----------------------------------
// Function: CTreeComps
//
// Returns: none
//
// Description: Default constructor for the gbm engine tree components.
//
// Parameters: double - the shrinkage parameter
//    unsigned long - number of training examples
//    unsigned long - number of features
//    double - fraction of training examples in bag
//    unsigned long - depth of the tree to grow
//    unsigned long - minimum number of observations in node
//    int - number of groups in data
//
//-----------------------------------
CTreeComps::CTreeComps(TreeParams treeConfig) :
  new_node_searcher_(treeConfig.depth, treeConfig.numberdatacolumns, treeConfig.min_obs_in_node),
  tree_(treeConfig.shrinkage, treeConfig.depth)
{
	this-> min_num_node_obs_ = treeConfig.min_obs_in_node;
	data_node_assignment_.resize(treeConfig.num_trainrows, 0);
}

//-----------------------------------
// Function: ~CTreeComps()
//
// Returns: none
//
// Description: default destructor for the gbmTreeComps
//
// Parameters: none
//
//-----------------------------------
CTreeComps::~CTreeComps()
{
}

//-----------------------------------
// Function: GrowTrees
//
// Returns: none
//
// Description: grows the tree
//
// Parameters: const CDataset ptr - pointer to the gbm data
//    int& - reference to  number of nodes in tree
//
//-----------------------------------
void CTreeComps::GrowTrees(const CDataset& data, double* adZ, const double* adFadj)
{
	#ifdef NOISY_DEBUG
	  Rprintf("Reset tree\n");
	#endif

	  //Reset tree and searcher
	  tree_.Reset();
	  new_node_searcher_.reset();

	#ifdef NOISY_DEBUG
	  Rprintf("grow tree\n");
	#endif

	tree_.grow(&(adZ[0]),
				data,
			&(adFadj[0]),
			min_num_node_obs_,
			data_node_assignment_,
			 new_node_searcher_);

	#ifdef NOISY_DEBUG
	  ptempTree->Print();
	#endif

	#ifdef NOISY_DEBUG
	  Rprintf("get node count=%d\n", tree_.GetNodeCount());
	#endif
}

//-----------------------------------
// Function: AdjustAndShrink
//
// Returns: none
//
// Description: adjusts the tree and shrinks.
//
// Parameters: none
//
//-----------------------------------
void CTreeComps::AdjustAndShrink(double * adFadj)
{
	tree_.Adjust(data_node_assignment_,
				  &(adFadj[0]),
				  min_num_node_obs_);

	#ifdef NOISY_DEBUG
	  tree_.Print();
	#endif
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
void CTreeComps::TransferTreeToRList(const CDataset &data,
	     int *aiSplitVar,
	     double *adSplitPoint,
	     int *aiLeftNode,
	     int *aiRightNode,
	     int *aiMissingNode,
	     double *adErrorReduction,
	     double *adWeight,
	     double *adPred,
	     VEC_VEC_CATEGORIES &vecSplitCodes,
	     int cCatSplitsOld)
{
	int nodeid = 0;
	if(tree_.GetRootNode())
	{
		tree_.GetRootNode()->TransferTreeToRList(nodeid,
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
											   tree_.GetShrinkageConst());
	}
	else
	{
	  throw GBM::Failure("Can't transfer to list - RootNode does not exist.");
	}
}

//-----------------------------------
// Function: PredictValid
//
// Returns: none
//
// Description: makes predictions on validation set
//
// Parameters: const CDataset ptr - ptr to gbm data
//
void CTreeComps::PredictValid(const CDataset& data, double* adFadj)
{
	tree_.PredictValid(data, data.get_validsize(), &(adFadj[0]));
}



