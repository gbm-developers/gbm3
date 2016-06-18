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
CNodeSearch::CNodeSearch(unsigned long treedepth, unsigned long minobs,
                         const parallel_details& parallel)
    : best_splits_(2 * treedepth + 1),
      num_terminal_nodes_(1),
      min_num_node_obs_(minobs),
      parallel_(parallel) {}

CNodeSearch::~CNodeSearch() {}

void CNodeSearch::GenerateAllSplits(vector<CNode*>& term_nodes_ptrs,
                                    const CDataset& kData, const Bag& kBag,
                                    const vector<double>& residuals,
                                    vector<unsigned long>& data_node_assigns) {
  const index_vector kColNumbers(kData.RandomOrder());
  VecNodeParams best_splits_updates(best_splits_);

#pragma omp parallel firstprivate(best_splits_updates) \
    num_threads(parallel_.get_num_threads())
  {
#pragma omp for schedule(static) nowait
    for (unsigned long ind = 0; ind < kData.get_num_features(); ++ind) {
      const int kVar = kColNumbers[ind];
      const int KVarClasses = kData.varclass(kVar);

      VecVarSplitters variable_splitters(num_terminal_nodes_, term_nodes_ptrs,
                                         min_num_node_obs_, ind, kVar,
                                         KVarClasses, kData.monotone(kVar));

      // distribute the observations in order to the correct node search
      for (unsigned long iOrderObs = 0; iOrderObs < kData.get_trainsize();
           iOrderObs++) {
        const unsigned long kWhichObs =
            kData.order_ptr()[kVar * kData.get_trainsize() + iOrderObs];
        if (kBag.get_element(kWhichObs)) {
          const int kNode = data_node_assigns[kWhichObs];
          const double kXVal = kData.x_value(kWhichObs, kVar);
          variable_splitters[kNode].IncorporateObs(
              kXVal, residuals[kWhichObs], kData.weight_ptr()[kWhichObs]);
        }
      }

      for (unsigned long node_num = 0; node_num < num_terminal_nodes_;
           node_num++) {
        variable_splitters[node_num].WrapUpCurrentVariable();
      }

      best_splits_updates += variable_splitters.proposal();
    }

#pragma omp critical
    { best_splits_ += best_splits_updates; }
  }
}

double CNodeSearch::CalcImprovementAndSplit(
    vector<CNode*>& term_nodes_ptrs, const CDataset& kData,
    vector<unsigned long>& data_node_assigns) {
  // search for the best split
  unsigned long bestnode = 0;
  double bestnode_improvement = -HUGE_VAL;
  for (unsigned long node_num = 0; node_num < num_terminal_nodes_; node_num++) {
    term_nodes_ptrs[node_num]->SetToSplit();
    if (best_splits_[node_num].get_improvement() > bestnode_improvement) {
      bestnode = node_num;
      bestnode_improvement = best_splits_[node_num].get_improvement();
    }
  }

  // Split Node if improvement is non-zero
  if (bestnode_improvement > 0.0) {
    // Split Node
    term_nodes_ptrs[bestnode]->SplitNode(best_splits_[bestnode]);
    num_terminal_nodes_ += 2;

    // Move kData to children nodes
    ReassignData(bestnode, term_nodes_ptrs, kData, data_node_assigns);

    // Add children to terminal node list
    term_nodes_ptrs[num_terminal_nodes_ - 2] =
        term_nodes_ptrs[bestnode]->right_child();
    term_nodes_ptrs[num_terminal_nodes_ - 1] =
        term_nodes_ptrs[bestnode]->missing_child();
    term_nodes_ptrs[bestnode] = term_nodes_ptrs[bestnode]->left_child();

    best_splits_[num_terminal_nodes_ - 2].ResetSplitProperties(
        term_nodes_ptrs[num_terminal_nodes_ - 2]->get_prediction() *
            term_nodes_ptrs[num_terminal_nodes_ - 2]->get_totalweight(),
        term_nodes_ptrs[num_terminal_nodes_ - 2]->get_totalweight(),
        term_nodes_ptrs[num_terminal_nodes_ - 2]->get_numobs());
    best_splits_[num_terminal_nodes_ - 1].ResetSplitProperties(
        term_nodes_ptrs[num_terminal_nodes_ - 1]->get_prediction() *
            term_nodes_ptrs[num_terminal_nodes_ - 1]->get_totalweight(),
        term_nodes_ptrs[num_terminal_nodes_ - 1]->get_totalweight(),
        term_nodes_ptrs[num_terminal_nodes_ - 1]->get_numobs());
    best_splits_[bestnode].ResetSplitProperties(
        term_nodes_ptrs[bestnode]->get_prediction() *
            term_nodes_ptrs[bestnode]->get_totalweight(),
        term_nodes_ptrs[bestnode]->get_totalweight(),
        term_nodes_ptrs[bestnode]->get_numobs());
  }

  return bestnode_improvement;
}

//----------------------------------------
// Function Members - Private
//----------------------------------------
void CNodeSearch::ReassignData(unsigned long splittednode_index,
                               vector<CNode*>& term_nodes_ptrs,
                               const CDataset& kData,
                               vector<unsigned long>& data_node_assigns) {
// assign observations to the correct node
#pragma omp parallel for schedule(static) \
    num_threads(parallel_.get_num_threads())
  for (unsigned long iObs = 0; iObs < kData.get_trainsize(); iObs++) {
    if (data_node_assigns[iObs] == splittednode_index) {
      signed char schWhichNode =
          term_nodes_ptrs[splittednode_index]->WhichNode(kData, iObs);
      if (schWhichNode == 1)  // goes right
      {
        data_node_assigns[iObs] = num_terminal_nodes_ - 2;
      } else if (schWhichNode == 0)  // is missing
      {
        data_node_assigns[iObs] = num_terminal_nodes_ - 1;
      }
      // those to the left stay with the same node assignment
    }
  }
}
