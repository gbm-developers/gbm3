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
CNodeSearch::CNodeSearch(unsigned long treedepth, unsigned long minobs):
variable_splitters_(2*treedepth+1, VarSplitter(minobs))
{
    num_terminal_nodes_ = 1;
    min_num_node_obs_ = minobs;
}


CNodeSearch::~CNodeSearch()
{
}

void CNodeSearch::GenerateAllSplits
(
		vector<CNode*>& term_nodes_ptrs,
		const CDataset& kData,
		double* residuals,
		vector<unsigned long>& data_node_assigns
)
{
	unsigned long kWhichObs = 0;
	const CDataset::index_vector kColNumbers(kData.RandomOrder());
	const CDataset::index_vector::const_iterator kFinalCol = kColNumbers.begin() + kData.get_num_features();

	for(CDataset::index_vector::const_iterator kIt=kColNumbers.begin();
	  kIt != kFinalCol;
	  kIt++)
	{
	  const int kVar = *kIt;
	  const int KVarClasses = kData.varclass(kVar);

	  for(unsigned long node_num=0; node_num < num_terminal_nodes_; node_num++)
	  {
		variable_splitters_[node_num].ResetForNewVar(kVar, KVarClasses);
	  }

	  // distribute the observations in order to the correct node search
	  for(unsigned long iOrderObs=0; iOrderObs < kData.get_trainsize(); iOrderObs++)
	  {
		  kWhichObs = kData.order_ptr()[kVar*kData.get_trainsize() + iOrderObs];
		  if(kData.get_bag_element(kWhichObs))
		  {
			  const int kNode = data_node_assigns[kWhichObs];
			  const double kXVal = kData.x_value(kWhichObs, kVar);
			  variable_splitters_[kNode].IncorporateObs(kXVal,
							residuals[kWhichObs],
							kData.weight_ptr()[kWhichObs],
							kData.monotone(kVar));
		  }
		}

		for(unsigned long node_num=0; node_num<num_terminal_nodes_; node_num++)
		{
			if(KVarClasses != 0) // evaluate if categorical split
			{
				variable_splitters_[node_num].EvaluateCategoricalSplit();
			}
				variable_splitters_[node_num].WrapUpCurrentVariable();
		}
	}

}


double CNodeSearch::CalcImprovementAndSplit
(
		vector<CNode*>& term_nodes_ptrs, const CDataset& kData,
		vector<unsigned long>& data_node_assigns
)
{
	// search for the best split
	unsigned long bestnode = 0;
	double bestnode_improvement = 0.0;
	for(unsigned long node_num=0; node_num < num_terminal_nodes_; node_num++)
	{
		variable_splitters_[node_num].SetToSplit();
		if(variable_splitters_[node_num].best_improvement() > bestnode_improvement)
		{
			bestnode = node_num;
			bestnode_improvement = variable_splitters_[node_num].best_improvement();
		}
	}


	// Split Node if improvement is non-zero
	if(bestnode_improvement != 0.0)
	{
		//Split Node
		variable_splitters_[bestnode].SetupNewNodes(*term_nodes_ptrs[bestnode]);
		num_terminal_nodes_ += 2;

		// Move kData to children nodes
		ReassignData(bestnode, term_nodes_ptrs, kData, data_node_assigns);

		// Add children to terminal node list
		term_nodes_ptrs[num_terminal_nodes_-2] = term_nodes_ptrs[bestnode]->right_node_ptr;
		term_nodes_ptrs[num_terminal_nodes_-1] = term_nodes_ptrs[bestnode]->missing_node_ptr;
		term_nodes_ptrs[bestnode] = term_nodes_ptrs[bestnode]->left_node_ptr;

		variable_splitters_[num_terminal_nodes_-2].Set(*term_nodes_ptrs[num_terminal_nodes_-2]);
		variable_splitters_[num_terminal_nodes_-1].Set(*term_nodes_ptrs[num_terminal_nodes_-1]);
		variable_splitters_[bestnode].Set(*term_nodes_ptrs[bestnode]);

	}

	return bestnode_improvement;
}


//----------------------------------------
// Function Members - Private
//----------------------------------------
void CNodeSearch::ReassignData
(
		unsigned long splittednode_index,
		vector<CNode*>& term_nodes_ptrs,
		const CDataset& kData,
		vector<unsigned long>& data_node_assigns
)
{
	// assign observations to the correct node
	for(unsigned long iObs=0; iObs < kData.get_trainsize(); iObs++)
	{
		if(data_node_assigns[iObs]==splittednode_index)
		{
		  signed char schWhichNode = term_nodes_ptrs[splittednode_index]->WhichNode(kData,iObs);
		  if(schWhichNode == 1) // goes right
		  {
			  data_node_assigns[iObs] = num_terminal_nodes_-2;
		  }
		  else if(schWhichNode == 0) // is missing
		  {
			  data_node_assigns[iObs] = num_terminal_nodes_-1;
		  }
		  // those to the left stay with the same node assignment
		}
	}
}

