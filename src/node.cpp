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
CNode::CNode(const NodeDef& kDefn) :
  prediction(kDefn.prediction()),
  totalweight(kDefn.get_totalweight()),
  numobs(kDefn.get_num_obs()),
  leftcategory() {

    splitvalue = 0.0;
    split_var = 0;
    improvement = 0.0;

    // Set children to NULL
	left_node_ptr = NULL;
	right_node_ptr = NULL;
	missing_node_ptr = NULL;

	// Set up split type and strategy
	splittype = kNone;
	node_strategy_ = new TerminalStrategy(this);

}

void CNode::SetStrategy()
{
	//delete nodeStrategy;
	switch(splittype)
	{
	case kNone:
		node_strategy_ = new TerminalStrategy(this);
		break;
	case kContinuous:
		node_strategy_ = new ContinuousStrategy(this);
		break;
	case kCategorical:
		node_strategy_ = new CategoricalStrategy(this);
		break;
	default:
		throw gbm_exception::Failure("Node State not recognised.");
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
    unsigned long min_num_node_obs
)
{
	node_strategy_->Adjust(min_num_node_obs);
}

void CNode::Predict
(
    const CDataset &kData,
    unsigned long rownum,
    double &delta_estimate
)
{
	node_strategy_->Predict(kData, rownum, delta_estimate);
}


void CNode::GetVarRelativeInfluence
(
    double* relative_influence
)
{
	node_strategy_->GetVarRelativeInfluence(relative_influence);
}

void CNode::PrintSubtree
(
 unsigned long indent
)
{
  node_strategy_->PrintSubTree(indent);
}

void CNode::SplitNode(NodeParams& childrenparams)
{

	// set up a continuous split
	if(childrenparams.split_class_==0)
	{
		splittype = kContinuous;
		SetStrategy();
	}
	else
	{
		splittype = kCategorical;
		SetStrategy();
		// the types are confused here
		leftcategory.resize(1 + (unsigned long)childrenparams.split_value_);
		std::copy(childrenparams.category_ordering_.begin(),
			  childrenparams.category_ordering_.begin() +
			  leftcategory.size(),
			  leftcategory.begin());
	}


	split_var = childrenparams.split_var_;
	splitvalue = childrenparams.split_value_;
	improvement = childrenparams.improvement_;

	left_node_ptr    = new CNode(childrenparams.left_);
	right_node_ptr   = new CNode(childrenparams.right_);
	missing_node_ptr = new CNode(childrenparams.missing_);



}

signed char CNode::WhichNode
(
    const CDataset &kData,
    unsigned long obs_num
)
{
	return node_strategy_->WhichNode(kData, obs_num);
}


void CNode::TransferTreeToRList
(
    int &node_id,
    const CDataset &kData,
    int* splivar,
    double* splitvalues,
    int* leftnodes,
    int* rightnodes,
    int* missingnodes,
    double* error_reduction,
    double* weights,
    double* predictions,
    VecOfVectorCategories &splitcodes_vec,
    int prev_categorical_split,
    double shrinkage
)
{
	node_strategy_->TransferTreeToRList(node_id,
										kData,
									splivar,
									splitvalues,
									leftnodes,
									rightnodes,
									missingnodes,
									error_reduction,
									weights,
									predictions,
									splitcodes_vec,
									prev_categorical_split,
									shrinkage);
}







