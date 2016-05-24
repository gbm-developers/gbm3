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
variable_splitters_(2*treeDepth+1, VarSplitter(minObs))
{
    num_terminal_nodes_ = 1;
    min_num_node_obs_ = minObs;
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
	const CDataset::index_vector colNumbers(data.RandomOrder());
	const CDataset::index_vector::const_iterator final = colNumbers.begin() + data.get_num_features();

	for(CDataset::index_vector::const_iterator it=colNumbers.begin();
	  it != final;
	  it++)
	{
	  const int iVar = *it;
	  const int cVarClasses = data.varclass(iVar);

	  for(unsigned long iNode=0; iNode < num_terminal_nodes_; iNode++)
	  {
		variable_splitters_[iNode].ResetForNewVar(iVar, cVarClasses);
	  }

	  // distribute the observations in order to the correct node search
	  for(unsigned long iOrderObs=0; iOrderObs < data.get_trainsize(); iOrderObs++)
	  {
		  iWhichObs = data.order_ptr()[iVar*data.get_trainsize() + iOrderObs];
		  if(data.get_bag_element(iWhichObs))
		  {
			  const int iNode = aiNodeAssign[iWhichObs];
			  const double dX = data.x_value(iWhichObs, iVar);
			  variable_splitters_[iNode].IncorporateObs(dX,
							adZ[iWhichObs],
							data.weight_ptr()[iWhichObs],
							data.monotone(iVar));
		  }
		}

		for(unsigned long iNode=0; iNode<num_terminal_nodes_; iNode++)
		{
			if(cVarClasses != 0) // evaluate if categorical split
			{
				variable_splitters_[iNode].EvaluateCategoricalSplit();
			}
				variable_splitters_[iNode].WrapUpCurrentVariable();
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
	for(unsigned long iNode=0; iNode < num_terminal_nodes_; iNode++)
	{
		variable_splitters_[iNode].SetToSplit();
		if(variable_splitters_[iNode].best_improvement() > dBestNodeImprovement)
		{
			iBestNode = iNode;
			dBestNodeImprovement = variable_splitters_[iNode].best_improvement();
		}
	}


	// Split Node if improvement is non-zero
	if(dBestNodeImprovement != 0.0)
	{
		//Split Node
		variable_splitters_[iBestNode].SetupNewNodes(*vecpTermNodes[iBestNode]);
		num_terminal_nodes_ += 2;

		// Move data to children nodes
		ReassignData(iBestNode, vecpTermNodes, data, aiNodeAssign);

		// Add children to terminal node list
		vecpTermNodes[num_terminal_nodes_-2] = vecpTermNodes[iBestNode]->right_node_ptr;
		vecpTermNodes[num_terminal_nodes_-1] = vecpTermNodes[iBestNode]->missing_node_ptr;
		vecpTermNodes[iBestNode] = vecpTermNodes[iBestNode]->left_node_ptr;

		variable_splitters_[num_terminal_nodes_-2].Set(*vecpTermNodes[num_terminal_nodes_-2]);
		variable_splitters_[num_terminal_nodes_-1].Set(*vecpTermNodes[num_terminal_nodes_-1]);
		variable_splitters_[iBestNode].Set(*vecpTermNodes[iBestNode]);

	}

	return dBestNodeImprovement;
}


//----------------------------------------
// Function Members - Private
//----------------------------------------
void CNodeSearch::ReassignData
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
			  aiNodeAssign[iObs] = num_terminal_nodes_-2;
		  }
		  else if(schWhichNode == 0) // is missing
		  {
			  aiNodeAssign[iObs] = num_terminal_nodes_-1;
		  }
		  // those to the left stay with the same node assignment
		}
	}
}

