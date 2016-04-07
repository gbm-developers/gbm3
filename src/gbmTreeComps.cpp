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
CTreeComps::CTreeComps(TreeParams treeConfig):
aNodeSearch(treeConfig.numColData, treeConfig.cMinObsInNode)
{
	this-> cMinObsInNode = treeConfig.cMinObsInNode;
	ptreeTemp = new CCARTTree(treeConfig.dShrinkage, treeConfig.cDepth);
	aiNodeAssign.resize(treeConfig.cTrain);

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
	delete ptreeTemp;
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
void CTreeComps::GrowTrees(const CDataset* pData, double* adZ, const double* adFadj)
{
	#ifdef NOISY_DEBUG
	  Rprintf("Reset tree\n");
	#endif

	  //Reset tree and searcher
	  ptreeTemp->Reset();
	  aNodeSearch.Reset();

	#ifdef NOISY_DEBUG
	  Rprintf("grow tree\n");
	#endif

	ptreeTemp->grow(&(adZ[0]),
	                *(pData),
	                &(adFadj[0]),
	                cMinObsInNode,
	                aiNodeAssign,
	                 aNodeSearch);

	#ifdef NOISY_DEBUG
	  ptempTree->Print();
	#endif

	#ifdef NOISY_DEBUG
	  Rprintf("get node count=%d\n", ptreeTemp->GetNodeCount());
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
	ptreeTemp->Adjust(aiNodeAssign,
	                  &(adFadj[0]),
	                  cMinObsInNode);

	#ifdef NOISY_DEBUG
	  ptreeTemp->Print();
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
void CTreeComps::TransferTreeToRList(const CDataset &pData,
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

	if(ptreeTemp->GetRootNode())
	{
		ptreeTemp->GetRootNode()->TransferTreeToRList(iNodeID,
													   pData,
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
													   ptreeTemp->GetShrinkageConst());
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
void CTreeComps::PredictValid(const CDataset* pData, double* adFadj)
{
	ptreeTemp->PredictValid(*(pData), pData->GetValidSize(), &(adFadj[0]));
}

//-----------------------------------
// Function: GetNodeAssign
//
// Returns: vector<unsigned long>
//
// Description: getter node assignments
//
// Parameters: none
//
//-----------------------------------
std::vector<unsigned long> CTreeComps::GetNodeAssign()
{
	return aiNodeAssign;
}

//-----------------------------------
// Function: GetTermNodes
//
// Returns: VEC_P_NODETERMINAL
//
// Description: getter for terminal nodes
//
// Parameters: none
//
//-----------------------------------
vector<CNode*> CTreeComps::GetTermNodes()
{
	return ptreeTemp->GetTermNodes();
}

//-----------------------------------
// Function: GetLambda
//
// Returns: double
//
// Description: get shrinkage
//
// Parameters: none
//
//-----------------------------------
const double CTreeComps::ShrinkageConstant() const
{
	return ptreeTemp->GetShrinkageConst();
}

//-----------------------------------
// Function: GetMinNodeObs
//
// Returns: unsigned long
//
// Description: get min no of observation in node
//
// Parameters: none
//
//-----------------------------------
unsigned long CTreeComps::GetMinNodeObs()
{
	return cMinObsInNode;
}


long CTreeComps::GetSizeOfTree()
{
	return ptreeTemp->GetNodeCount();
}
const long CTreeComps::GetSizeOfTree() const
{
	return ptreeTemp->GetNodeCount();
}


