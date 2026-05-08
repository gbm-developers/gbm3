//-----------------------------------
//
// File: adaboost.cpp
//
// Description: distribution used for adaboosting.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "adaboost.h"
#include <algorithm>
#include <cmath>
#include <memory>

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

double exp_diff(double log_x, double log_y) {
  if (log_x == log_y) {
    return 0.0;
  }

  const double max_log = std::max(log_x, log_y);
  return std::exp(max_log) *
         (std::exp(log_x - max_log) - std::exp(log_y - max_log));
}
}

//----------------------------------------
// Function Members - Private
//----------------------------------------
CAdaBoost::CAdaBoost(const parallel_details& parallel)
    : CDistribution(parallel) {}

//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CAdaBoost::Create(DataDistParams& distparams) {
  return new CAdaBoost(distparams.parallel);
}

CAdaBoost::~CAdaBoost() {}

void CAdaBoost::ComputeWorkingResponse(const CDataset& kData, const Bag& kBag,
                                       const double* kFuncEstimate,
                                       std::vector<double>& residuals) {
#pragma omp parallel for schedule(static, get_array_chunk_size()) \
  num_threads(get_num_threads())
  for (unsigned long i = 0; i < kData.get_trainsize(); i++) {
    residuals[i] = -(2 * kData.y_ptr()[i] - 1) *
                   std::exp(-(2 * kData.y_ptr()[i] - 1) *
                            (kData.offset_ptr()[i] + kFuncEstimate[i]));
  }
}

double CAdaBoost::InitF(const CDataset& kData) {
  double log_numerator = -HUGE_VAL;
  double log_denominator = -HUGE_VAL;

#pragma omp parallel num_threads(get_num_threads())
  {
    double thread_log_numerator = -HUGE_VAL;
    double thread_log_denominator = -HUGE_VAL;

#pragma omp for schedule(static, get_array_chunk_size())
    for (unsigned long i = 0; i < kData.get_trainsize(); i++) {
      if (kData.weight_ptr()[i] > 0.0) {
        const double log_weight = std::log(kData.weight_ptr()[i]);
        if (kData.y_ptr()[i] == 1.0) {
          thread_log_numerator =
              log_add_exp(thread_log_numerator,
                          log_weight - kData.offset_ptr()[i]);
        } else {
          thread_log_denominator =
              log_add_exp(thread_log_denominator,
                          log_weight + kData.offset_ptr()[i]);
        }
      }
    }

#pragma omp critical
    {
      log_numerator = log_add_exp(log_numerator, thread_log_numerator);
      log_denominator = log_add_exp(log_denominator, thread_log_denominator);
    }
  }

  return 0.5 * (log_numerator - log_denominator);
}

double CAdaBoost::Deviance(const CDataset& kData, const Bag& kBag,
                           const double* kFuncEstimate) {
  double loss = 0.0;
  double weight = 0.0;

  // Switch to validation set if necessary
  unsigned long num_of_rows_in_set = kData.get_size_of_set();

#pragma omp parallel for schedule(static, get_array_chunk_size()) \
    reduction(+ : loss, weight) num_threads(get_num_threads())
  for (unsigned long i = 0; i < num_of_rows_in_set; i++) {
    loss += kData.weight_ptr()[i] *
            std::exp(-(2 * kData.y_ptr()[i] - 1) *
                     (kData.offset_ptr()[i] + kFuncEstimate[i]));
    weight += kData.weight_ptr()[i];
  }

  // TODO: Check if weights are all zero for validation set
  if ((weight == 0.0) && (loss == 0.0)) {
    return nan("");
  } else if (weight == 0.0) {
    return HUGE_VAL;
  }

  return loss / weight;
}

void CAdaBoost::FitBestConstant(const CDataset& kData, const Bag& kBag,
                                const double* kFuncEstimate,
                                unsigned long num_terminalnodes,
                                std::vector<double>& residuals,
                                CCARTTree& tree) {
  unsigned long obs_num = 0;
  unsigned long node_num = 0;
  numerator_bestconstant_.resize(num_terminalnodes);
  numerator_bestconstant_.assign(numerator_bestconstant_.size(), 0.0);
  denominator_bestconstant_.resize(num_terminalnodes);
  denominator_bestconstant_.assign(denominator_bestconstant_.size(), 0.0);
  std::vector<double> scale_bestconstant_(num_terminalnodes, -HUGE_VAL);

  for (obs_num = 0; obs_num < kData.get_trainsize(); obs_num++) {
    if (kBag.get_element(obs_num) && kData.weight_ptr()[obs_num] > 0.0) {
      const unsigned long node = tree.get_node_assignments()[obs_num];
      const double outcome = 2 * kData.y_ptr()[obs_num] - 1;
      const double deltafunc_est =
          kFuncEstimate[obs_num] + kData.offset_ptr()[obs_num];
      const double log_weight =
          std::log(kData.weight_ptr()[obs_num]) - outcome * deltafunc_est;

      double scaled_weight = 1.0;
      if (log_weight > scale_bestconstant_[node]) {
        if (scale_bestconstant_[node] != -HUGE_VAL) {
          const double rescale = std::exp(scale_bestconstant_[node] - log_weight);
          numerator_bestconstant_[node] *= rescale;
          denominator_bestconstant_[node] *= rescale;
        }
        scale_bestconstant_[node] = log_weight;
      } else {
        scaled_weight = std::exp(log_weight - scale_bestconstant_[node]);
      }

      numerator_bestconstant_[node] += outcome * scaled_weight;
      denominator_bestconstant_[node] += scaled_weight;
    }
  }

  for (node_num = 0; node_num < num_terminalnodes; node_num++) {
    if (tree.has_node(node_num)) {
      if (denominator_bestconstant_[node_num] == 0) {
        tree.get_terminal_nodes()[node_num]->set_prediction(0.0);
      } else {
        tree.get_terminal_nodes()[node_num]->set_prediction(
            numerator_bestconstant_[node_num] /
            denominator_bestconstant_[node_num]);
      }
    }
  }
}

double CAdaBoost::BagImprovement(const CDataset& kData, const Bag& kBag,
                                 const double* kFuncEstimate,
                                 const double kShrinkage,
                                 const std::vector<double>& kDeltaEstimate) {
  double returnvalue = 0.0;
  double weight = 0.0;

#pragma omp parallel for schedule(static, get_array_chunk_size()) \
    reduction(+ : returnvalue, weight) num_threads(get_num_threads())
  for (unsigned long i = 0; i < kData.get_trainsize(); i++) {
    if (!kBag.get_element(i)) {
      const double func_est = kFuncEstimate[i] + kData.offset_ptr()[i];
      const double outcome = 2 * kData.y_ptr()[i] - 1;
      const double log_loss = -outcome * func_est;
      const double log_adjusted_loss =
          -outcome * (func_est + kShrinkage * kDeltaEstimate[i]);

      returnvalue +=
          kData.weight_ptr()[i] *
          exp_diff(log_loss, log_adjusted_loss);
      weight += kData.weight_ptr()[i];
    }
  }

  return returnvalue / weight;
}
