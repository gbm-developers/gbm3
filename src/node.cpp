//  GBM by Greg Ridgeway  Copyright (C) 2003

//-----------------------------------
// Includes
//-----------------------------------
#include "node.h"
#include "terminal_strategy.h"
#include "continuous_strategy.h"
#include "categorical_strategy.h"

//----------------------------------------
// Function Members - Public
//----------------------------------------
CNode::CNode(const NodeDef& defn) :
  prediction(defn.prediction()),
  totalweight(defn.get_totalweight()),
  numobs(defn.get_num_obs()),
  leftcategory() {

    splitvalue = 0.0;
    split_var = 0;
    improvement = 0.0;

    // Set children to NULL
	left_node_ptr = NULL;
	right_node_ptr = NULL;
	missing_node_ptr = NULL;

	// Set up split type and strategy
	splittype = none;
	node_strategy_ = new TerminalStrategy(this);

}

void CNode::SetStrategy()
{
	//delete nodeStrategy;
	switch(splittype)
	{
	case none:
		node_strategy_ = new TerminalStrategy(this);
		break;
	case continuous:
		node_strategy_ = new ContinuousStrategy(this);
		break;
	case categorical:
		node_strategy_ = new CategoricalStrategy(this);
		break;
	default:
		throw GBM::Failure("Node State not recognised.");
		break;
	}
}

CNode::~CNode()
{
	// Each node is responsible for deleting its
	// children and its strategy
    delete left_node_ptr;
    delete right_node_ptr;
    delete missing_node_ptr;
    delete node_strategy_;
}

void CNode::Adjust
(
    unsigned long cMinObsInNode
)
{
	node_strategy_->Adjust(cMinObsInNode);
}

void CNode::Predict
(
    const CDataset &data,
    unsigned long iRow,
    double &dFadj
)
{
	node_strategy_->Predict(data, iRow, dFadj);
}


void CNode::GetVarRelativeInfluence
(
    double *adRelInf
)
{
	node_strategy_->GetVarRelativeInfluence(adRelInf);
}

void CNode::PrintSubtree
(
 unsigned long cIndent
)
{
  node_strategy_->PrintSubTree(cIndent);
}

void CNode::SplitNode(NodeParams& childrenParams)
{

	// set up a continuous split
	if(childrenParams.split_class_==0)
	{
		splittype = continuous;
		SetStrategy();
	}
	else
	{
		splittype = categorical;
		SetStrategy();
		// the types are confused here
		leftcategory.resize(1 + (ULONG)childrenParams.split_value_);
		std::copy(childrenParams.category_ordering_.begin(),
			  childrenParams.category_ordering_.begin() +
			  leftcategory.size(),
			  leftcategory.begin());
	}


	split_var = childrenParams.split_var_;
	splitvalue = childrenParams.split_value_;
	improvement = childrenParams.improvement_;

	left_node_ptr    = new CNode(childrenParams.left_);
	right_node_ptr   = new CNode(childrenParams.right_);
	missing_node_ptr = new CNode(childrenParams.missing_);



}

signed char CNode::WhichNode
(
    const CDataset &data,
    unsigned long iObs
)
{
	return node_strategy_->WhichNode(data, iObs);
}


void CNode::TransferTreeToRList
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
	node_strategy_->TransferTreeToRList(iNodeID,
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







