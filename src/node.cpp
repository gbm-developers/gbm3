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
CNode::CNode(double nodePrediction,
		double trainingWeight, long numObs):aiLeftCategory()
{
    dPrediction = nodePrediction;
    dTrainW = trainingWeight;
    cN = numObs;

    dSplitValue = 0.0;
    iSplitVar = 0;
    dImprovement = 0.0;

    // Set children to NULL
	pLeftNode = NULL;
	pRightNode = NULL;
	pMissingNode = NULL;

	// Set up split type and strategy
	splitType = none;
	nodeStrategy = new TerminalStrategy(this);

}

void CNode::SetStrategy()
{
	delete nodeStrategy;
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

void CNode::SplitNode()
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
		aiLeftCategory.resize(1 + (ULONG)childrenParams.SplitValue);
					  std::copy(childrenParams.aiBestCategory.begin(),
								childrenParams.aiBestCategory.begin() + aiLeftCategory.size(),
								 aiLeftCategory.begin());
	}

	iSplitVar = childrenParams.SplitVar;
	dSplitValue = childrenParams.SplitValue;
	dImprovement = childrenParams.ImprovedResiduals;

	pLeftNode    = new CNode(childrenParams.LeftWeightResiduals/childrenParams.LeftTotalWeight, childrenParams.LeftTotalWeight,
									childrenParams.LeftNumObs);
	pRightNode   = new CNode(childrenParams.RightWeightResiduals/childrenParams.RightTotalWeight,
							childrenParams.RightTotalWeight, childrenParams.RightNumObs);
	pMissingNode = new CNode(childrenParams.MissingWeightResiduals/childrenParams.MissingTotalWeight,
							childrenParams.MissingTotalWeight, childrenParams.MissingNumObs);

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







