//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//
//  File:       node_search.cpp
//
//------------------------------------------------------------------------------
//-----------------------------------
// Includes
//-----------------------------------
#include "node_search.h"

//----------------------------------------
// Function Members - Public
//----------------------------------------
CNodeSearch::CNodeSearch(int treeDepth, int numColData, unsigned long minObs):
variableSplitters(2*treeDepth+1, VarSplitter(minObs))//variableSplitters(numColData, VarSplitter(minObs))
{
    cTerminalNodes = 1;
    minNumObs = minObs;
    totalCache = 0;

}


CNodeSearch::~CNodeSearch()
{
}

void CNodeSearch::GenerateAllSplits
(
		vector<CNode*>& vecpTermNodes,
		const CDataset& data,
		double* adZ,
		vector<unsigned long>& aiNodeAssign
)
{
	unsigned long iWhichObs = 0;
	const CDataset::index_vector colNumbers(data.random_order());
	const CDataset::index_vector::const_iterator final = colNumbers.begin() + data.get_numFeatures();

	for(long iNode = 0; iNode < cTerminalNodes; iNode++)
	{
		// If node has split then skip
		if(vecpTermNodes[iNode]->splitAssigned) continue;
		variableSplitters[iNode].Set(*vecpTermNodes[iNode]);

	}

	  for(CDataset::index_vector::const_iterator it=colNumbers.begin();
	      it != final;
	      it++)
	    {
	      const int iVar = *it;
	      const int cVarClasses = data.varclass(iVar);

	      for(long iNode=0; iNode < cTerminalNodes; iNode++)
	        {
	    	  if(vecpTermNodes[iNode]->splitAssigned) continue;
	    	  variableSplitters[iNode].ResetForNewVar(iVar, cVarClasses);
	        }

	      // distribute the observations in order to the correct node search
	      for(long iOrderObs=0; iOrderObs < data.get_trainSize(); iOrderObs++)
	        {
		  iWhichObs = data.order_ptr()[iVar*data.get_trainSize() + iOrderObs];

		  if(vecpTermNodes[aiNodeAssign[iWhichObs]]->splitAssigned) continue;
		  if(data.GetBagElem(iWhichObs))
	            {
		      const int iNode = aiNodeAssign[iWhichObs];
		      const double dX = data.x_value(iWhichObs, iVar);
		      variableSplitters[iNode].IncorporateObs(dX,
							adZ[iWhichObs],
							data.weight_ptr()[iWhichObs],
							data.monotone(iVar));
	            }
	        }
	        for(long iNode=0; iNode<cTerminalNodes; iNode++)
	        {
	        	if(vecpTermNodes[iNode]->splitAssigned) continue;
	            if(cVarClasses != 0) // evaluate if categorical split
	            {
		      variableSplitters[iNode].EvaluateCategoricalSplit();
	            }
	            variableSplitters[iNode].WrapUpCurrentVariable();
	        }
	    }
}


double CNodeSearch::CalcImprovementAndSplit
(
		vector<CNode*>& vecpTermNodes, const CDataset& data,
		vector<unsigned long>& aiNodeAssign
)
{
	// search for the best split
	long iBestNode = 0;
	double dBestNodeImprovement = 0.0;
	for(long iNode=0; iNode < cTerminalNodes; iNode++)
	{

		if(variableSplitters[iNode].BestImprovement() > dBestNodeImprovement)
		{
			iBestNode = iNode;
			dBestNodeImprovement = variableSplitters[iNode].BestImprovement();
			vecpTermNodes[iNode]->childrenParams =  variableSplitters[iNode].GetBestSplit();
		}

		if(vecpTermNodes[iNode]->splitAssigned) continue;
		vecpTermNodes[iNode]->SplitAssign();

	}
	// Split Node if improvement is non-zero
	if(dBestNodeImprovement != 0.0)
	{
		//Split Node
		//vecpTermNodes[iBestNode]->SplitNode();
		variableSplitters[iBestNode].SetupNewNodes(*vecpTermNodes[iBestNode]);
		cTerminalNodes += 2;

		// Move data to children nodes
		ReAssignData(iBestNode, vecpTermNodes, data, aiNodeAssign);

		// Add children to terminal node list
		vecpTermNodes[cTerminalNodes-2] = vecpTermNodes[iBestNode]->pRightNode;
		vecpTermNodes[cTerminalNodes-1] = vecpTermNodes[iBestNode]->pMissingNode;
		vecpTermNodes[iBestNode] = vecpTermNodes[iBestNode]->pLeftNode;

	}

	return dBestNodeImprovement;
}


//----------------------------------------
// Function Members - Private
//----------------------------------------
void CNodeSearch::ReAssignData
(
		long splittedNodeIndex,
		vector<CNode*>& vecpTermNodes,
		const CDataset& data,
		vector<unsigned long>& aiNodeAssign
)
{
	// assign observations to the correct node
	for(long iObs=0; iObs < data.get_trainSize(); iObs++)
	{
		if(aiNodeAssign[iObs]==splittedNodeIndex)
		{
		  signed char schWhichNode = vecpTermNodes[splittedNodeIndex]->WhichNode(data,iObs);
		  if(schWhichNode == 1) // goes right
		  {
			  aiNodeAssign[iObs] = cTerminalNodes-2;
		  }
		  else if(schWhichNode == 0) // is missing
		  {
			  aiNodeAssign[iObs] = cTerminalNodes-1;
		  }
		  // those to the left stay with the same node assignment
		  }
	}
}



