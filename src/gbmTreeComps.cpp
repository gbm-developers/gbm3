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
#include "gbmTreeComps.h"

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
  aNodeSearch(treeConfig.cDepth, treeConfig.numColData, treeConfig.cMinObsInNode),
  tree(treeConfig.dShrinkage, treeConfig.cDepth)
{
	this-> cMinObsInNode = treeConfig.cMinObsInNode;
	aiNodeAssign.resize(treeConfig.cTrain, 0);
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
	  tree.Reset();
	  aNodeSearch.reset();

	#ifdef NOISY_DEBUG
	  Rprintf("grow tree\n");
	#endif

	tree.grow(&(adZ[0]),
				data,
			&(adFadj[0]),
			cMinObsInNode,
			aiNodeAssign,
			 aNodeSearch);

	#ifdef NOISY_DEBUG
	  ptempTree->Print();
	#endif

	#ifdef NOISY_DEBUG
	  Rprintf("get node count=%d\n", tree.GetNodeCount());
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
	tree.Adjust(aiNodeAssign,
				  &(adFadj[0]),
				  cMinObsInNode);

	#ifdef NOISY_DEBUG
	  tree.Print();
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
	int iNodeID = 0;

	if(tree.GetRootNode())
	{
		tree.GetRootNode()->TransferTreeToRList(iNodeID,
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
											   tree.GetShrinkageConst());
	}
	else
	{
	  throw GBM::failure("Can't transfer to list - RootNode does not exist.");
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
	tree.PredictValid(data, data.get_validsize(), &(adFadj[0]));
}



