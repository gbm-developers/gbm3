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
depthOfTree(depth), shrinkageConst(shrinkage)
{
    pRootNode = NULL;
    cTotalNodeCount = 1;
    dError = 0.0;

    // Calculate original
}


CCARTTree::~CCARTTree()
{
	delete pRootNode;
}

void CCARTTree::Reset()
{
  delete pRootNode;
  pRootNode = NULL;
  vecpTermNodes.resize(2*depthOfTree + 1, NULL);
  cTotalNodeCount = 1;
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
	 (depthOfTree < 1))
	{
	  throw GBM::invalid_argument();
	}

  double dSumZ = 0.0;
  double dSumZ2 = 0.0;
  double dTotalW = 0.0;

#ifdef NOISY_DEBUG
  Rprintf("initial tree calcs\n");
#endif

  	  // Move to data -- FOR TIME BEING
	for(long iObs=0; iObs<data.get_trainSize(); iObs++)
	{
		// aiNodeAssign tracks to which node each training obs belongs
		aiNodeAssign[iObs] = 0;

		if(data.GetBagElem(iObs))
		{
			// get the initial sums and sum of squares and total weight
			dSumZ += data.weight_ptr()[iObs]*adZ[iObs];
			dSumZ2 += data.weight_ptr()[iObs]*adZ[iObs]*adZ[iObs];
			dTotalW += data.weight_ptr()[iObs];
		}
	}

  dError = dSumZ2-dSumZ*dSumZ/dTotalW;
  pRootNode = new CNode(NodeDef(dSumZ, dTotalW, data.GetTotalInBag()));
  vecpTermNodes[0] = pRootNode;
  aNodeSearch.SetRootNode(*pRootNode);

  // build the tree structure
#ifdef NOISY_DEBUG
  Rprintf("Building tree 1 ");
#endif

  for(long cDepth=0; cDepth < depthOfTree; cDepth++)
  {
#ifdef NOISY_DEBUG
      Rprintf("%d ",cDepth);
#endif
      
    // Generate all splits
    aNodeSearch.GenerateAllSplits(vecpTermNodes, data, &(adZ[0]), aiNodeAssign);
    double bestImprov = aNodeSearch.CalcImprovementAndSplit(vecpTermNodes, data, aiNodeAssign);

    // Make the best split if possible
	if(bestImprov == 0.0)
	{
	  break;
	}
	// setup the new nodes and add them to the tree

	cTotalNodeCount += 3;

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
  int i=0;
  
  for(i=data.nrow() - nValid; i<data.nrow(); i++)
    {
      pRootNode->Predict(data, i, adFadj[i]);
      adFadj[i] *= shrinkageConst;
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

	pRootNode->Adjust(cMinObsInNode);

	// predict for the training observations
	for(iObs=0; iObs<aiNodeAssign.size(); iObs++)
	{
		adFadj[iObs] = vecpTermNodes[aiNodeAssign[iObs]]->dPrediction;
	}
}


void CCARTTree::Print()
{
    if(pRootNode)
    {
      pRootNode->PrintSubtree(0);
      Rprintf("shrinkage: %f\n",shrinkageConst);
      Rprintf("initial error: %f\n\n",dError);
    }
}


CNode* CCARTTree::GetRootNode()
{
	return pRootNode;
}

const CNode* CCARTTree::GetRootNode() const
{
	return pRootNode;
}


