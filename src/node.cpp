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
  dTrainW(defn.totalWeight),
  cN(defn.numObs),
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
	/*signed char schWhichNode;
	switch(splitType)
	{
	case none:
		dFadj = dPrediction;
		break;
	case continuous:
		schWhichNode = WhichNode(data,iRow);
		if(schWhichNode == -1)
		{
		  pLeftNode->Predict(data, iRow, dFadj);
		}
		else if(schWhichNode == 1)
		{
		  pRightNode->Predict(data, iRow, dFadj);
		}
		else
		{
		  pMissingNode->Predict(data, iRow, dFadj);
		}
		break;
	case categorical:
		schWhichNode = WhichNode(data,iRow);
		if(schWhichNode == -1)
		{
		  pLeftNode->Predict(data, iRow, dFadj);
		}
		else if(schWhichNode == 1)
		{
		  pRightNode->Predict(data, iRow, dFadj);
		}
		else
		{
		  pMissingNode->Predict(data, iRow, dFadj);
		}
		break;
	default:
			throw GBM::failure("Node State not recognised.");
			break;
	}*/
	nodeStrategy->Predict(data, iRow, dFadj);
}


void CNode::GetVarRelativeInfluence
(
    double *adRelInf
)
{
	/*switch(splitType)
	{
	case none:
		return;
		break;
	case continuous:
		adRelInf[iSplitVar] += dImprovement;
		pLeftNode->GetVarRelativeInfluence(adRelInf);
		pRightNode->GetVarRelativeInfluence(adRelInf);
		break;
	case categorical:
		adRelInf[iSplitVar] += dImprovement;
		pLeftNode->GetVarRelativeInfluence(adRelInf);
		pRightNode->GetVarRelativeInfluence(adRelInf);
		break;
	default:
			throw GBM::failure("Node State not recognised.");
			break;
	}*/
	nodeStrategy->GetVarRelativeInfluence(adRelInf);
}

void CNode::PrintSubtree
(
 unsigned long cIndent
)
{
	/*const std::size_t cLeftCategory = aiLeftCategory.size();
	switch(splitType)
	{
	case none:
		for(long i=0; i< cIndent; i++) Rprintf("  ");
				  Rprintf("N=%f, Prediction=%f *\n",
					 dTrainW,
					 dPrediction);
		break;
	case continuous:

		for(long i=0; i< cIndent; i++) Rprintf("  ");
		Rprintf("N=%f, Improvement=%f, Prediction=%f, NA pred=%f\n",
		  dTrainW,
		 dImprovement,
		  dPrediction,
		  (pMissingNode == NULL ? 0.0 : pMissingNode->dPrediction));

		for(long i=0; i< cIndent; i++) Rprintf("  ");
		Rprintf("V%d in ", iSplitVar);
		for(long i=0; i<cLeftCategory; i++)
		  {
			Rprintf("%d", aiLeftCategory[i]);
			if(i<cLeftCategory-1) Rprintf(",");
		  }
		Rprintf("\n");
		pLeftNode->PrintSubtree((cIndent+1));

		for(long i=0; i< cIndent; i++) Rprintf("  ");
		Rprintf("V%d not in ", iSplitVar);
		for(long i=0; i<cLeftCategory; i++)
		  {
			Rprintf("%d", aiLeftCategory[i]);
			if(i<cLeftCategory-1) Rprintf(",");
		  }
		Rprintf("\n");
		pRightNode->PrintSubtree(cIndent+1);

		for(long i=0; i< cIndent; i++) Rprintf("  ");
		Rprintf("missing\n");
		pMissingNode->PrintSubtree(cIndent+1);
		break;
	case categorical:

		for(long i=0; i< cIndent; i++) Rprintf("  ");
		Rprintf("N=%f, Improvement=%f, Prediction=%f, NA pred=%f\n",
		  dTrainW,
		  dImprovement,
		  dPrediction,
		  (pMissingNode == NULL ? 0.0 : pMissingNode->dPrediction));

		for(long i=0; i< cIndent; i++) Rprintf("  ");
		Rprintf("V%d in ",iSplitVar);
		for(long i=0; i<cLeftCategory; i++)
		  {
			Rprintf("%d", aiLeftCategory[i]);
			if(i<cLeftCategory-1) Rprintf(",");
		  }
		Rprintf("\n");
		pLeftNode->PrintSubtree((cIndent+1));

		for(long i=0; i< cIndent; i++) Rprintf("  ");
		Rprintf("V%d not in ", iSplitVar);
		for(long i=0; i<cLeftCategory; i++)
		  {
			Rprintf("%d", aiLeftCategory[i]);
			if(i<cLeftCategory-1) Rprintf(",");
		  }
		Rprintf("\n");
		pRightNode->PrintSubtree(cIndent+1);

		for(long i=0; i< cIndent; i++) Rprintf("  ");
		Rprintf("missing\n");
		pMissingNode->PrintSubtree(cIndent+1);
		break;
	default:
			throw GBM::failure("Node State not recognised.");
			break;
	}*/
  nodeStrategy->PrintSubTree(cIndent);
}

void CNode::SplitAssign()
{
	splitAssigned = true;
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
	/*signed char ReturnValue = 0;
	double dX = data.x_value(iObs, iSplitVar);
	switch(splitType)
	{
	case none:
		 if(!ISNA(dX))
			{
				if(dX < dSplitValue)
				{
					ReturnValue = -1;
				}
				else
				{
					ReturnValue = 1;
				}
			}
			// if missing value returns 0


		break;
	case continuous:
		 if(!ISNA(dX))
			{
				if(dX < dSplitValue)
				{
					ReturnValue = -1;
				}
				else
				{
					ReturnValue = 1;
				}
			}
			// if missing value returns 0

			return ReturnValue;
		break;
	case categorical:
		if(!ISNA(dX))
		{
		  if(std::find(aiLeftCategory.begin(),
			   aiLeftCategory.end(),
			   (ULONG)dX) != aiLeftCategory.end())
			{
				ReturnValue = -1;
			}
			else
			{
				ReturnValue = 1;
			}
		}
		// if missing value returns 0
		break;
	default:
			throw GBM::failure("Node State not recognised.");
			break;
	}

	return ReturnValue;*/
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
	/*switch(splitType)
	{
	case none:
	{
		int iThisNodeID = iNodeID;
		aiSplitVar[iNodeID] = -1;
		adSplitPoint[iNodeID] = dShrinkage*dPrediction;
		aiLeftNode[iNodeID] = -1;
		aiRightNode[iNodeID] = -1;
		aiMissingNode[iNodeID] = -1;
		adErrorReduction[iNodeID] = 0.0;
		adWeight[iNodeID] = dTrainW;
		adPred[iNodeID] = dShrinkage*dPrediction;

		iNodeID++;
		break;
	}
	case continuous:
	{
		int iThisNodeID = iNodeID;
		aiSplitVar[iThisNodeID] = iSplitVar;
		adSplitPoint[iThisNodeID] = dSplitValue;
		adErrorReduction[iThisNodeID] = dImprovement;
		adWeight[iThisNodeID] = dTrainW;
		adPred[iThisNodeID] = dShrinkage*dPrediction;

		iNodeID++;
		aiLeftNode[iThisNodeID] = iNodeID;
		pLeftNode->TransferTreeToRList(iNodeID,
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

		aiRightNode[iThisNodeID] = iNodeID;
		pRightNode->TransferTreeToRList(iNodeID,
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
		aiMissingNode[iThisNodeID] = iNodeID;
		pMissingNode->TransferTreeToRList(iNodeID,
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
		break;
	}
	case categorical:
	{
		int iThisNodeID = iNodeID;
		unsigned long cCatSplits = vecSplitCodes.size();
		unsigned long i = 0;
		int cLevels = data.varclass(iSplitVar);
		const std::size_t cLeftCategory = aiLeftCategory.size();

		aiSplitVar[iThisNodeID] = iSplitVar;
		adSplitPoint[iThisNodeID] = cCatSplits+cCatSplitsOld; // 0 based
		adErrorReduction[iThisNodeID] = dImprovement;
		adWeight[iThisNodeID] = dTrainW;
		adPred[iThisNodeID] = dShrinkage*dPrediction;

		vecSplitCodes.push_back(VEC_CATEGORIES());

		vecSplitCodes[cCatSplits].resize(cLevels,1);
		for(i=0; i<cLeftCategory; i++)
		  {
			vecSplitCodes[cCatSplits][aiLeftCategory[i]] = -1;
		  }

		iNodeID++;
		aiLeftNode[iThisNodeID] = iNodeID;
		pLeftNode->TransferTreeToRList(iNodeID,
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
		aiRightNode[iThisNodeID] = iNodeID;
		pRightNode->TransferTreeToRList(iNodeID,
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

		aiMissingNode[iThisNodeID] = iNodeID;
		pMissingNode->TransferTreeToRList(iNodeID,
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
		break;
	}
	default:
	{
		throw GBM::failure("Node State not recognised.");
		break;
	}
	}*/
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







