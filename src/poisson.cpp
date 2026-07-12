//-----------------------------------
//
// File: poisson.cpp
//
// Description: poisson distribution for GBM.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "poisson.h"
#include <cmath>

namespace {
double log_add_exp(double total, double value) {
  if (total == -HUGE_VAL) {
    return value;
  }
  if (value == -HUGE_VAL) {
    return total;
  }
  if (value > total) {
    return value + std::log1p(std::exp(total - value));
  }
  return total + std::log1p(std::exp(value - total));
}
}

//----------------------------------------
// Function Members - Private
//----------------------------------------
CPoisson::CPoisson(const parallel_details& parallel)
    : CDistribution(parallel) {}

//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CPoisson::Create(DataDistParams& distparams) {
  return new CPoisson(distparams.parallel);
}

CPoisson::~CPoisson() {}

void CPoisson::ComputeWorkingResponse(const CDataset& kData, const Bag& kBag,
                                      const double* kFuncEstimate,
                                      std::vector<double>& residuals) {
// compute working response
#pragma omp parallel for schedule(static, get_array_chunk_size()) \
  num_threads(get_num_threads())
  for (unsigned long i = 0; i < kData.get_trainsize(); i++) {
    const double delta_func_est = kFuncEstimate[i] + kData.offset_ptr()[i];
    residuals[i] = kData.y_ptr()[i] - std::exp(delta_func_est);
  }
}

double CPoisson::InitF(const CDataset& kData) {
  double log_numerator = -HUGE_VAL;
  double log_denominator = -HUGE_VAL;
  double min = -19.0;
  double max = 19.0;
  double max_offset = -HUGE_VAL;
  double min_offset = HUGE_VAL;

  for (unsigned long i = 0; i < kData.get_trainsize(); i++) {
    if (kData.weight_ptr()[i] > 0.0) {
      const double log_weight = std::log(kData.weight_ptr()[i]);
      log_denominator =
          log_add_exp(log_denominator, log_weight + kData.offset_ptr()[i]);
      if (kData.y_ptr()[i] > 0.0) {
        log_numerator =
            log_add_exp(log_numerator, log_weight + std::log(kData.y_ptr()[i]));
      }
      max_offset = R::fmax2(max_offset, kData.offset_ptr()[i]);
      min_offset = R::fmin2(min_offset, kData.offset_ptr()[i]);
    }
  }

  if (log_numerator == -HUGE_VAL) {
    return min;
  }

  double init_func_est = log_numerator - log_denominator;

  // Keep eta_i = offset_i + init_func_est within [min, max] for every
  // observation, not just the offset-free init_func_est itself.
  if (max_offset + init_func_est > max) {
    init_func_est = max - max_offset;
  }
  if (min_offset + init_func_est < min) {
    init_func_est = min - min_offset;
  }

  return init_func_est;
}

double CPoisson::Deviance(const CDataset& kData, const Bag& kBag,
                          const double* kFuncEstimate) {
  double loss = 0.0;
  double weight = 0.0;

  // Switch to validation set if necessary
  unsigned long num_rows_in_set = kData.get_size_of_set();

#pragma omp parallel for schedule(static, get_array_chunk_size()) \
    reduction(+ : loss, weight) num_threads(get_num_threads())
  for (unsigned long i = 0; i < num_rows_in_set; i++) {
    loss += kData.weight_ptr()[i] *
            (kData.y_ptr()[i] * (kData.offset_ptr()[i] + kFuncEstimate[i]) -
             std::exp(kData.offset_ptr()[i] + kFuncEstimate[i]));
    weight += kData.weight_ptr()[i];
  }

  // TODO: Check if weights are all zero for validation set
  if ((weight == 0.0) && (loss == 0.0)) {
    return nan("");
  } else if (weight == 0.0) {
    return copysign(HUGE_VAL, -loss);
  }
  return -2 * loss / weight;
}

void CPoisson::FitBestConstant(const CDataset& kData, const Bag& kBag,
                               const double* kFuncEstimate,
                               unsigned long num_terminalnodes,
                               std::vector<double>& residuals,
                               CCARTTree& tree) {
  unsigned long obs_num = 0;
  unsigned long node_num = 0;
  std::vector<double> numerator_vec(num_terminalnodes, 0.0);
  std::vector<double> denominator_vec(num_terminalnodes, 0.0);
  std::vector<double> max_vec(num_terminalnodes, -HUGE_VAL);
  std::vector<double> min_vec(num_terminalnodes, HUGE_VAL);

  for (obs_num = 0; obs_num < kData.get_trainsize(); obs_num++) {
    if (kBag.get_element(obs_num)) {
      const double delta_func_est =
          kData.offset_ptr()[obs_num] + kFuncEstimate[obs_num];
      const unsigned long node_num = tree.get_node_assignments()[obs_num];
      numerator_vec[tree.get_node_assignments()[obs_num]] +=
          kData.weight_ptr()[obs_num] * kData.y_ptr()[obs_num];
      denominator_vec[tree.get_node_assignments()[obs_num]] +=
          kData.weight_ptr()[obs_num] * std::exp(delta_func_est);

      max_vec[node_num] = R::fmax2(delta_func_est, max_vec[node_num]);
      min_vec[node_num] = R::fmin2(delta_func_est, min_vec[node_num]);
    }
  }

  for (node_num = 0; node_num < num_terminalnodes; node_num++) {
    if (tree.has_node(node_num)) {
      if (numerator_vec[node_num] == 0.0) {
        // DEBUG: if vecdNum==0 then prediction = -Inf
        // Not sure what else to do except plug in an arbitrary
        //   negative number, -1? -10? Let's use -1, then make
        //   sure |adF| < 19 always.
        tree.get_terminal_nodes()[node_num]->set_prediction(-19.0);
      } else if (denominator_vec[node_num] == 0.0) {
        tree.get_terminal_nodes()[node_num]->set_prediction(0.0);
      } else {
        tree.get_terminal_nodes()[node_num]->set_prediction(
            std::log(numerator_vec[node_num] / denominator_vec[node_num]));
      }
      tree.get_terminal_nodes()[node_num]->set_prediction(
          R::fmin2(tree.get_terminal_nodes()[node_num]->get_prediction(),
                   19 - max_vec[node_num]));
      tree.get_terminal_nodes()[node_num]->set_prediction(
          R::fmax2(tree.get_terminal_nodes()[node_num]->get_prediction(),
                   -19 - min_vec[node_num]));
    }
  }
}

double CPoisson::BagImprovement(const CDataset& kData, const Bag& kBag,
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
          kData.weight_ptr()[i] *
          (kData.y_ptr()[i] * kShrinkage * kDeltaEstimate[i] -
           std::exp(delta_func_est + kShrinkage * kDeltaEstimate[i]) +
           std::exp(delta_func_est));
      weight += kData.weight_ptr()[i];
    }
  }

  return returnvalue / weight;
}
