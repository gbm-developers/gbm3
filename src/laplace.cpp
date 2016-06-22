//-----------------------------------
//
// File: laplace.cpp
//
// Description: laplace distribution for GBM.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "laplace.h"
#include <vector>

//----------------------------------------
// Function Members - Private
//----------------------------------------
CLaplace::CLaplace(const parallel_details& parallel)
    : CDistribution(parallel), mpLocM_("Other") {}

//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CLaplace::Create(DataDistParams& distparams) {
  return new CLaplace(distparams.parallel);
}

CLaplace::~CLaplace() {}

void CLaplace::ComputeWorkingResponse(const CDataset& kData, const Bag& kBag,
                                      const double* kFuncEstimate,
                                      std::vector<double>& residuals) {
#pragma omp parallel for schedule(static, get_array_chunk_size()) \
  num_threads(get_num_threads())
  for (unsigned long i = 0; i < kData.get_trainsize(); i++) {
    residuals[i] =
        (kData.y_ptr()[i] - kData.offset_ptr()[i] - kFuncEstimate[i]) > 0.0
            ? 1.0
            : -1.0;
  }
}

double CLaplace::InitF(const CDataset& kData) {
  std::vector<double> arr(kData.get_trainsize());

#pragma omp parallel for schedule(static, get_array_chunk_size()) \
  num_threads(get_num_threads())
  for (unsigned long ii = 0; ii < kData.get_trainsize(); ii++) {
    arr[ii] = kData.y_ptr()[ii] - kData.offset_ptr()[ii];
  }

  return mpLocM_.WeightedQuantile(kData.get_trainsize(), &arr[0],
                                  kData.weight_ptr(), 0.5);  // median
}

double CLaplace::Deviance(const CDataset& kData, const Bag& kBag,
                          const double* kFuncEstimates) {
  double loss = 0.0;
  double weight = 0.0;

  unsigned long num_rows_in_set = kData.get_size_of_set();

#pragma omp parallel for schedule(static, get_array_chunk_size()) \
    reduction(+ : loss, weight) num_threads(get_num_threads())
  for (unsigned long i = 0; i < num_rows_in_set; i++) {
    loss += kData.weight_ptr()[i] *
            fabs(kData.y_ptr()[i] - kData.offset_ptr()[i] - kFuncEstimates[i]);
    weight += kData.weight_ptr()[i];
  }

  // TODO: Check if weights are all zero for validation set
  if ((weight == 0.0) && (loss == 0.0)) {
    return nan("");
  } else if (weight == 0.0) {
    return copysign(HUGE_VAL, loss);
  }

  return loss / weight;
}

// DEBUG: needs weighted median
void CLaplace::FitBestConstant(const CDataset& kData, const Bag& kBag,
                               const double* kFuncEstimate,
                               unsigned long num_terminal_nodes,
                               std::vector<double>& residuals,
                               CCARTTree& tree) {
  unsigned long node_num = 0;
  unsigned long obs_num = 0;
  unsigned long vec_num = 0;
  double offset;

  std::vector<double> adArr(kData.get_trainsize());
  std::vector<double> adW2(kData.get_trainsize());

  for (node_num = 0; node_num < num_terminal_nodes; node_num++) {
    if (tree.get_terminal_nodes()[node_num]->get_numobs() >=
        tree.min_num_obs_required()) {
      vec_num = 0;
      for (obs_num = 0; obs_num < kData.get_trainsize(); obs_num++) {
        if (kBag.get_element(obs_num) &&
            (tree.get_node_assignments()[obs_num] == node_num)) {
          offset = kData.offset_ptr()[obs_num];
          adArr[vec_num] =
              kData.y_ptr()[obs_num] - offset - kFuncEstimate[obs_num];
          adW2[vec_num] = kData.weight_ptr()[obs_num];
          vec_num++;
        }
      }

      tree.get_terminal_nodes()[node_num]->set_prediction(
          mpLocM_.WeightedQuantile(vec_num, &adArr[0], &adW2[0],
                                   0.5));  // median
    }
  }
}

double CLaplace::BagImprovement(const CDataset& kData, const Bag& kBag,
                                const double* kFuncEstimate,
                                const double kShrinkage,
                                const std::vector<double>& kDeltaEstimate) {
  double returnvalue = 0.0;
  double weight = 0.0;

#pragma omp parallel for schedule(static, get_array_chunk_size()) \
    reduction(+ : returnvalue, weight) num_threads(get_num_threads())
  for (unsigned long i = 0; i < kData.get_trainsize(); i++) {
    if (!kBag.get_element(i)) {
      const double delta_func_est = kFuncEstimate[i] + kData.offset_ptr()[i];

      returnvalue +=
          kData.weight_ptr()[i] * (fabs(kData.y_ptr()[i] - delta_func_est) -
                                   fabs(kData.y_ptr()[i] - delta_func_est -
                                        kShrinkage * kDeltaEstimate[i]));
      weight += kData.weight_ptr()[i];
    }
  }

  return returnvalue / weight;
}
