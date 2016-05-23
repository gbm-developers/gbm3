//  GBM by Greg Ridgeway  Copyright (C) 2003

//-----------------------------------
// Includes
//-----------------------------------
#include "node.h"
#include "terminalStrategy.h"
#include "continuousStrategy.h"
#include "categoricalStrategy.h"

//----------------------------------------
// Function Members - Public
//----------------------------------------
CNode::CNode(const NodeDef& defn) :
  dPrediction(defn.prediction()),
  dTrainW(defn.get_totalweight()),
  cN(defn.get_num_obs()),
  aiLeftCategory() {

    dSplitValue = 0.0;
    iSplitVar = 0;
    dImprovement = 0.0;

    // Set children to NULL
	pLeftNode = NULL;
	pRightNode = NULL;
	pMissingNode = NULL;

	// Set up split type and strategy
	splitType = none;
	splitAssigned = false;
	nodeStrategy = new TerminalStrategy(this);

}

void CNode::SetStrategy()
{
	//delete nodeStrategy;
	switch(splitType)
	{
	case none:
		nodeStrategy = new TerminalStrategy(this);
		break;
	case continuous:
		nodeStrategy = new ContinuousStrategy(this);
		break;
	case categorical:
		nodeStrategy = new CategoricalStrategy(this);
		break;
	default:
		throw GBM::failure("Node State not recognised.");
		break;
	}
}

CNode::~CNode()
{
	// Each node is responsible for deleting its
	// children and its strategy
    delete pLeftNode;
    delete pRightNode;
    delete pMissingNode;
    delete nodeStrategy;
}

void CNode::Adjust
(
    unsigned long cMinObsInNode
)
{
	/*switch(splitType)
	{
	case none:
		return;
		break;
	case continuous:
		pLeftNode->Adjust(cMinObsInNode);
		pRightNode->Adjust(cMinObsInNode);

		if((pMissingNode->splitType == none) && (pMissingNode->cN < cMinObsInNode))
		{
			dPrediction = ((pLeftNode->dTrainW)*(pLeftNode->dPrediction) +
				 (pRightNode->dTrainW)*(pRightNode->dPrediction))/
			(pLeftNode->dTrainW + pRightNode->dTrainW);
			pMissingNode->dPrediction = dPrediction;
		}
		else
		{
			pMissingNode->Adjust(cMinObsInNode);
			dPrediction =
			((pLeftNode->dTrainW)*(pLeftNode->dPrediction) +
			(pRightNode->dTrainW)*  (pRightNode->dPrediction) +
			(pMissingNode->dTrainW)*(pMissingNode->dPrediction))/
			(pLeftNode->dTrainW + pRightNode->dTrainW + pMissingNode->dTrainW);
		}
		break;
	case categorical:
		pLeftNode->Adjust(cMinObsInNode);
		pRightNode->Adjust(cMinObsInNode);

		if((pMissingNode->splitType == none) && (pMissingNode->cN < cMinObsInNode))
		{
			dPrediction = ((pLeftNode->dTrainW)*(pLeftNode->dPrediction) +
				 (pRightNode->dTrainW)*(pRightNode->dPrediction))/
			(pLeftNode->dTrainW + pRightNode->dTrainW);
			pMissingNode->dPrediction = dPrediction;
		}
		else
		{
			pMissingNode->Adjust(cMinObsInNode);
			dPrediction =
			((pLeftNode->dTrainW)*(pLeftNode->dPrediction) +
			(pRightNode->dTrainW)*  (pRightNode->dPrediction) +
			(pMissingNode->dTrainW)*(pMissingNode->dPrediction))/
			(pLeftNode->dTrainW + pRightNode->dTrainW + pMissingNode->dTrainW);
		}
		break;
	default:
			throw GBM::failure("Node State not recognised.");
			break;
	}*/
	nodeStrategy->Adjust(cMinObsInNode);
}

void CNode::Predict
(
    const CDataset &data,
    unsigned long iRow,
    double &dFadj
)
{
	nodeStrategy->Predict(data, iRow, dFadj);
}


void CNode::GetVarRelativeInfluence
(
    double *adRelInf
)
{
	nodeStrategy->GetVarRelativeInfluence(adRelInf);
}

void CNode::PrintSubtree
(
 unsigned long cIndent
)
{
  nodeStrategy->PrintSubTree(cIndent);
}

void CNode::SplitNode(NodeParams& childrenParams)
{

	// set up a continuous split
	if(childrenParams.SplitClass==0)
	{
		splitType = continuous;
		SetStrategy();
	}
	else
	{
		splitType = categorical;
		SetStrategy();
		// the types are confused here
		aiLeftCategory.resize(1 + (ULONG)childrenParams.SplitValue);
		std::copy(childrenParams.aiBestCategory.begin(),
			  childrenParams.aiBestCategory.begin() +
			  aiLeftCategory.size(),
			  aiLeftCategory.begin());
	}


	iSplitVar = childrenParams.SplitVar;
	dSplitValue = childrenParams.SplitValue;
	dImprovement = childrenParams.ImprovedResiduals;

	pLeftNode    = new CNode(childrenParams.left);
	pRightNode   = new CNode(childrenParams.right);
	pMissingNode = new CNode(childrenParams.missing);



}

signed char CNode::WhichNode
(
    const CDataset &data,
    unsigned long iObs
)
{
	return nodeStrategy->WhichNode(data, iObs);
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
	nodeStrategy->TransferTreeToRList(iNodeID,
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







