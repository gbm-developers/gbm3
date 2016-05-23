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
variableSplitters(2*treeDepth+1, VarSplitter(minObs))
{
    cTerminalNodes = 1;
    minNumObs = minObs;
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
	const CDataset::index_vector::const_iterator final = colNumbers.begin() + data.get_num_features();

	for(CDataset::index_vector::const_iterator it=colNumbers.begin();
	  it != final;
	  it++)
	{
	  const int iVar = *it;
	  const int cVarClasses = data.varclass(iVar);

	  for(unsigned long iNode=0; iNode < cTerminalNodes; iNode++)
	  {
		variableSplitters[iNode].ResetForNewVar(iVar, cVarClasses);
	  }

	  // distribute the observations in order to the correct node search
	  for(unsigned long iOrderObs=0; iOrderObs < data.get_trainsize(); iOrderObs++)
	  {
		  iWhichObs = data.order_ptr()[iVar*data.get_trainsize() + iOrderObs];
		  if(data.get_bag_element(iWhichObs))
		  {
			  const int iNode = aiNodeAssign[iWhichObs];
			  const double dX = data.x_value(iWhichObs, iVar);
			  variableSplitters[iNode].IncorporateObs(dX,
							adZ[iWhichObs],
							data.weight_ptr()[iWhichObs],
							data.monotone(iVar));
		  }
		}

		for(unsigned long iNode=0; iNode<cTerminalNodes; iNode++)
		{
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
	unsigned long iBestNode = 0;
	double dBestNodeImprovement = 0.0;
	for(unsigned long iNode=0; iNode < cTerminalNodes; iNode++)
	{
		variableSplitters[iNode].SetToSplit();
		if(variableSplitters[iNode].best_improvement() > dBestNodeImprovement)
		{
			iBestNode = iNode;
			dBestNodeImprovement = variableSplitters[iNode].best_improvement();
		}
	}


	// Split Node if improvement is non-zero
	if(dBestNodeImprovement != 0.0)
	{
		//Split Node
		variableSplitters[iBestNode].SetupNewNodes(*vecpTermNodes[iBestNode]);
		cTerminalNodes += 2;

		// Move data to children nodes
		ReAssignData(iBestNode, vecpTermNodes, data, aiNodeAssign);

		// Add children to terminal node list
		vecpTermNodes[cTerminalNodes-2] = vecpTermNodes[iBestNode]->pRightNode;
		vecpTermNodes[cTerminalNodes-1] = vecpTermNodes[iBestNode]->pMissingNode;
		vecpTermNodes[iBestNode] = vecpTermNodes[iBestNode]->pLeftNode;

		variableSplitters[cTerminalNodes-2].Set(*vecpTermNodes[cTerminalNodes-2]);
		variableSplitters[cTerminalNodes-1].Set(*vecpTermNodes[cTerminalNodes-1]);
		variableSplitters[iBestNode].Set(*vecpTermNodes[iBestNode]);

	}

	return dBestNodeImprovement;
}


//----------------------------------------
// Function Members - Private
//----------------------------------------
void CNodeSearch::ReAssignData
(
		unsigned long splittedNodeIndex,
		vector<CNode*>& vecpTermNodes,
		const CDataset& data,
		vector<unsigned long>& aiNodeAssign
)
{
	// assign observations to the correct node
	for(unsigned long iObs=0; iObs < data.get_trainsize(); iObs++)
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

