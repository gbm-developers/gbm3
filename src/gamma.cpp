//-----------------------------------
//
// File: gamma.cpp
//
// Description: Gamma distribution with natural
//  log link function ( mean = std::exp(prediction) )
//
// Notes: -19 <= prediction <= +19
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------

#include "gamma.h"
#include <cmath>

//----------------------------------------
// Function Members - Private
//----------------------------------------
CGamma::CGamma(const parallel_details& parallel) : CDistribution(parallel) {}

//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CGamma::Create(DataDistParams& distparams) {
  return new CGamma(distparams.parallel);
}

CGamma::~CGamma() {}

void CGamma::ComputeWorkingResponse(const CDataset& kData, const Bag& kBag,
                                    const double* kFuncEstimate,
                                    std::vector<double>& residuals) {
  if (!(kData.y_ptr() && kFuncEstimate && &(residuals[0]) &&
        kData.weight_ptr())) {
    throw gbm_exception::InvalidArgument();
  }

#pragma omp parallel for schedule(static) num_threads(get_num_threads())
  for (unsigned long i = 0; i < kData.get_trainsize(); i++) {
    const double deltafunc_est = kFuncEstimate[i] + kData.offset_ptr()[i];
    residuals[i] = kData.y_ptr()[i] * std::exp(-deltafunc_est) - 1.0;
  }
}

double CGamma::InitF(const CDataset& kData) {
  double sum = 0.0;
  double totalweight = 0.0;
  double min = -19.0;
  double max = +19.0;
  double initfunc_est = 0.0;

#pragma omp parallel for schedule(static) \
    reduction(+ : sum, totalweight) num_threads(get_num_threads())
  for (unsigned long i = 0; i < kData.get_trainsize(); i++) {
    sum += kData.weight_ptr()[i] * kData.y_ptr()[i] *
           std::exp(-kData.offset_ptr()[i]);
    totalweight += kData.weight_ptr()[i];
  }

  if (sum <= 0.0) {
    initfunc_est = min;
  } else {
    initfunc_est = std::log(sum / totalweight);
  }

  if (initfunc_est < min) {
    initfunc_est = min;
  }
  if (initfunc_est > max) {
    initfunc_est = max;
  }
  return initfunc_est;
}

double CGamma::Deviance(const CDataset& kData, const Bag& kBag,
                        const double* kFuncEstimate) {
  double loss = 0.0;
  double weight = 0.0;

  unsigned long num_rows_in_set = kData.get_size_of_set();

#pragma omp parallel for schedule(static) \
    reduction(+ : loss, weight) num_threads(get_num_threads())
  for (unsigned long i = 0; i < num_rows_in_set; i++) {
    const double deltafunc_est = kFuncEstimate[i] + kData.offset_ptr()[i];
    loss += kData.weight_ptr()[i] *
            (kData.y_ptr()[i] * std::exp(-deltafunc_est) + deltafunc_est);
    weight += kData.weight_ptr()[i];
  }

  // TODO: Check if weights are all zero for validation set
  if ((weight == 0.0) && (loss == 0.0)) {
    return nan("");
  } else if (weight == 0.0) {
    return copysign(HUGE_VAL, loss);
  }

  return 2 * loss / weight;
}

void CGamma::FitBestConstant(const CDataset& kData, const Bag& kBag,
                             const double* kFuncEstimate,
                             unsigned long num_terminalnodes,
                             std::vector<double>& residuals, CCARTTree& tree) {
  double deltafunc_estimate = 0.0;
  unsigned long obs_num = 0;
  unsigned long node_num = 0;
  double maxval = 19.0;
  double minval = -19.0;

  vector<double> numerator_vec(num_terminalnodes, 0.0);
  vector<double> denominator_vec(num_terminalnodes, 0.0);
  vector<double> max_vec(num_terminalnodes, -HUGE_VAL);
  vector<double> min_vec(num_terminalnodes, HUGE_VAL);

  for (obs_num = 0; obs_num < kData.get_trainsize(); obs_num++) {
    if (kBag.get_element(obs_num)) {
      deltafunc_estimate = kFuncEstimate[obs_num] + kData.offset_ptr()[obs_num];
      numerator_vec[tree.get_node_assignments()[obs_num]] +=
          kData.weight_ptr()[obs_num] * kData.y_ptr()[obs_num] *
          std::exp(-deltafunc_estimate);
      denominator_vec[tree.get_node_assignments()[obs_num]] +=
          kData.weight_ptr()[obs_num];

      // Keep track of largest and smallest prediction in each node
      max_vec[tree.get_node_assignments()[obs_num]] = R::fmax2(
          deltafunc_estimate, max_vec[tree.get_node_assignments()[obs_num]]);
      min_vec[tree.get_node_assignments()[obs_num]] = R::fmin2(
          deltafunc_estimate, min_vec[tree.get_node_assignments()[obs_num]]);
    }
  }

  for (node_num = 0; node_num < num_terminalnodes; node_num++) {
    if (tree.get_terminal_nodes()[node_num] != NULL) {
      if (numerator_vec[node_num] == 0.0) {
        // Taken from poisson.cpp

        // DEBUG: if numerator_vec==0 then prediction = -Inf
        // Not sure what else to do except plug in an arbitrary
        //   negative number, -1? -10? Let's use -19, then make
        //   sure |kFuncEstimate| < 19 always.
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

double CGamma::BagImprovement(const CDataset& kData, const Bag& kBag,
                              const double* kFuncEstimate,
                              const double kShrinkage,
                              const std::vector<double>& kDeltaEstimate) {
  double returnvalue = 0.0;
  double weight = 0.0;

#pragma omp parallel for schedule(static) \
    reduction(+ : returnvalue, weight) num_threads(get_num_threads())
  for (unsigned long i = 0; i < kData.get_trainsize(); i++) {
    if (!kBag.get_element(i)) {
      const double deltafunc_est = kFuncEstimate[i] + kData.offset_ptr()[i];
      returnvalue += kData.weight_ptr()[i] *
                     (kData.y_ptr()[i] * std::exp(-deltafunc_est) *
                          (1.0 - exp(-kShrinkage * kDeltaEstimate[i])) -
                      kShrinkage * kDeltaEstimate[i]);
      weight += kData.weight_ptr()[i];
    }
  }

  return 2 * returnvalue / weight;
}
