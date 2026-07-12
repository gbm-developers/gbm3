//-----------------------------------
//
// File: tweedie.cpp
//
// Description: Tweedie distribution with natural
//  log link function ( mean = std::exp(prediction) )
//
// Notes: -19 <= prediction <= +19, Parameter dPower defaults to 1.5
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "tweedie.h"
#include <cmath>
#include <typeinfo>
#include <iostream>
#include <vector>
#include <deque>
#include <fstream>

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
CTweedie::CTweedie(double power) { power_ = power; }

//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CTweedie::Create(DataDistParams& distparams) {
  // Extract misc from second column of response
  double power = Rcpp::as<double>(distparams.misc[0]);
  if (!gbm_functions::has_value(power)) {
    throw gbm_exception::Failure(
        "Tweedie distribution requires misc to initialization.");
  }
  return new CTweedie(power);
}

CTweedie::~CTweedie() {}

void CTweedie::ComputeWorkingResponse(const CDataset& kData, const Bag& kBag,
                                      const double* kFuncEstimates,
                                      std::vector<double>& residuals) {
  unsigned long i = 0;
  double delta_func_est = 0.0;

  if (!(kData.y_ptr() && kFuncEstimates &&
        kData.weight_ptr())) {
    throw gbm_exception::InvalidArgument();
  }

  for (i = 0; i < kData.get_trainsize(); i++) {
    delta_func_est = kFuncEstimates[i] + kData.offset_ptr()[i];
    residuals[i] =
        kData.y_ptr()[i] * std::exp(delta_func_est * (1.0 - power_)) -
        exp(delta_func_est * (2.0 - power_));
  }
}

double CTweedie::InitF(const CDataset& kData) {
  double log_numerator = -HUGE_VAL;
  double log_denominator = -HUGE_VAL;
  double min = -19.0;
  double max = +19.0;
  unsigned long i = 0;
  double init_func_est = 0.0;
  double max_offset = -HUGE_VAL;
  double min_offset = HUGE_VAL;

  for (i = 0; i < kData.get_trainsize(); i++) {
    if (kData.weight_ptr()[i] > 0.0) {
      const double log_weight = std::log(kData.weight_ptr()[i]);
      log_denominator =
          log_add_exp(log_denominator,
                      log_weight + kData.offset_ptr()[i] * (2.0 - power_));
      if (kData.y_ptr()[i] > 0.0) {
        log_numerator =
            log_add_exp(log_numerator,
                        log_weight + std::log(kData.y_ptr()[i]) +
                            kData.offset_ptr()[i] * (1.0 - power_));
      }
      max_offset = R::fmax2(max_offset, kData.offset_ptr()[i]);
      min_offset = R::fmin2(min_offset, kData.offset_ptr()[i]);
    }
  }

  if (log_numerator == -HUGE_VAL) {
    init_func_est = min;
  } else {
    init_func_est = log_numerator - log_denominator;
  }

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

double CTweedie::Deviance(const CDataset& kData, const Bag& kBag,
                          const double* kFuncEstimate) {
  double delta_func_est = 0.0;
  unsigned long i = 0;
  double loss = 0.0;
  double weight = 0.0;

  // Switch to validation set if necessary
  unsigned long num_rows_in_set = kData.get_size_of_set();

  for (i = 0; i < num_rows_in_set; i++) {
    delta_func_est = kFuncEstimate[i] + kData.offset_ptr()[i];
    loss += kData.weight_ptr()[i] *
            (pow(kData.y_ptr()[i], 2.0 - power_) /
                 ((1.0 - power_) * (2.0 - power_)) -
             kData.y_ptr()[i] * std::exp(delta_func_est * (1.0 - power_)) /
                 (1.0 - power_) +
             exp(delta_func_est * (2.0 - power_)) / (2.0 - power_));
    weight += kData.weight_ptr()[i];
  }

  // TODO: Check if weights are all zero for validation set
  if ((weight == 0.0) && (loss == 0.0)) {
    return nan("");
  } else if (weight == 0.0) {
    return copysign(HUGE_VAL, loss);
  }

  return 2.0 * loss / weight;
}

void CTweedie::FitBestConstant(const CDataset& kData, const Bag& kBag,
                               const double* kFuncEstimate,
                               unsigned long num_terminalnodes,
                               std::vector<double>& residuals,
                               CCARTTree& tree) {
  double delta_func_est = 0.0;
  unsigned long obs_num = 0;
  unsigned long node_num = 0;
  double maxval = 19.0;
  double minval = -19.0;

  std::vector<double> numerator_vec(num_terminalnodes, 0.0);
  std::vector<double> denominator_vec(num_terminalnodes, 0.0);
  std::vector<double> max_vec(num_terminalnodes, -HUGE_VAL);
  std::vector<double> min_vec(num_terminalnodes, HUGE_VAL);

  for (obs_num = 0; obs_num < kData.get_trainsize(); obs_num++) {
    if (kBag.get_element(obs_num)) {
      delta_func_est = kFuncEstimate[obs_num] + kData.offset_ptr()[obs_num];
      numerator_vec[tree.get_node_assignments()[obs_num]] +=
          kData.weight_ptr()[obs_num] * kData.y_ptr()[obs_num] *
          std::exp(delta_func_est * (1.0 - power_));
      denominator_vec[tree.get_node_assignments()[obs_num]] +=
          kData.weight_ptr()[obs_num] *
          std::exp(delta_func_est * (2.0 - power_));

      // Keep track of largest and smallest prediction in each node
      max_vec[tree.get_node_assignments()[obs_num]] = R::fmax2(
          delta_func_est, max_vec[tree.get_node_assignments()[obs_num]]);
      min_vec[tree.get_node_assignments()[obs_num]] = R::fmin2(
          delta_func_est, min_vec[tree.get_node_assignments()[obs_num]]);
    }
  }

  for (node_num = 0; node_num < num_terminalnodes; node_num++) {
    if (tree.has_node(node_num)) {
      if (numerator_vec[node_num] == 0.0) {
        // Taken from poisson.cpp

        // DEBUG: if vecdNum==0 then prediction = -Inf
        // Not sure what else to do except plug in an arbitrary
        //   negative number, -1? -10? Let's use -19, then make
        //   sure |adF| < 19 always.
        tree.get_terminal_nodes()[node_num]->set_prediction(minval);
      }

      else if (denominator_vec[node_num] == 0.0) {
        tree.get_terminal_nodes()[node_num]->set_prediction(0.0);
      }

      else {
        tree.get_terminal_nodes()[node_num]->set_prediction(
            std::log(numerator_vec[node_num] / denominator_vec[node_num]));
      }

      if (max_vec[node_num] +
              tree.get_terminal_nodes()[node_num]->get_prediction() >
          maxval) {
        tree.get_terminal_nodes()[node_num]->set_prediction(maxval -
                                                            max_vec[node_num]);
      }
      if (min_vec[node_num] +
              tree.get_terminal_nodes()[node_num]->get_prediction() <
          minval) {
        tree.get_terminal_nodes()[node_num]->set_prediction(minval -
                                                            min_vec[node_num]);
      }
    }
  }
}

double CTweedie::BagImprovement(const CDataset& kData, const Bag& kBag,
                                const double* kFuncEstimate,
                                const double kShrinkage,
                                const std::vector<double>& kDeltaEstimate) {
  double returnvalue = 0.0;
  double delta_func_estimate = 0.0;
  double weight = 0.0;
  unsigned long i = 0;

  for (i = 0; i < kData.get_trainsize(); i++) {
    if (!kBag.get_element(i)) {
      delta_func_estimate = kFuncEstimate[i] + kData.offset_ptr()[i];

      returnvalue +=
          kData.weight_ptr()[i] *
          (std::exp(delta_func_estimate * (1.0 - power_)) * kData.y_ptr()[i] /
               (1.0 - power_) *
               (std::exp(kShrinkage * kDeltaEstimate[i] * (1.0 - power_)) -
                1.0) +
           std::exp(delta_func_estimate * (2.0 - power_)) / (2.0 - power_) *
               (1.0 - exp(kShrinkage * kDeltaEstimate[i] * (2.0 - power_))));
      weight += kData.weight_ptr()[i];
    }
  }

  return 2.0 * returnvalue / weight;
}
