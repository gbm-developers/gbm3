//-----------------------------------
//
// File: bernoulli.cpp
//
// Description: bernoulli distribution.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------

#include "bernoulli.h"
#include <memory>

//----------------------------------------
// Function Members - Private
//----------------------------------------
CBernoulli::CBernoulli()
    : terminalnode_capped_(false), terminalnode_cap_level_(10) {}

//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CBernoulli::Create(DataDistParams& distparams) {
  return new CBernoulli();
}

CBernoulli::~CBernoulli() {}

void CBernoulli::ComputeWorkingResponse(const CDataset& kData, const Bag& kBag,
                                        const double* kFuncEstimate,
                                        std::vector<double>& residuals) {
  double prob = 0.0;
  double deltafunc_est = 0.0;

  for (unsigned long i = 0; i < kData.get_trainsize(); i++) {
    deltafunc_est = kFuncEstimate[i] + kData.offset_ptr()[i];
    prob = 1.0 / (1.0 + std::exp(-deltafunc_est));

    residuals[i] = kData.y_ptr()[i] - prob;
  }
}

double CBernoulli::InitF(const CDataset& kData) {
  // Newton method for solving for F
  // should take about 3-6 iterations.

  double initfunc_est = 0.0;
  double newtonstep = 1.0;  // set to 1 initially

  for (int itcount = 0; (itcount < 6) && (std::abs(newtonstep) > 0.001);
       ++itcount) {
    double numerator = 0.0;
    double denominator = 0.0;

    for (unsigned long i = 0; i < kData.get_trainsize(); i++) {
      const double dTemp =
          1.0 / (1.0 + std::exp(-(kData.offset_ptr()[i] + initfunc_est)));
      numerator += kData.weight_ptr()[i] * (kData.y_ptr()[i] - dTemp);
      denominator += kData.weight_ptr()[i] * dTemp * (1.0 - dTemp);
    }
    newtonstep = numerator / denominator;
    initfunc_est += newtonstep;
  }

  return initfunc_est;
}

double CBernoulli::Deviance(const CDataset& kData, const Bag& kBag,
                            const double* kFuncEstimate) {
  unsigned long i = 0;
  double loss = 0.0;
  double deltafunc_est = 0.0;
  double weight = 0.0;

  // Switch to validation set if necessary
  unsigned long num_of_rows_in_set = kData.get_size_of_set();

  for (i = 0; i != num_of_rows_in_set; i++) {
    deltafunc_est = kFuncEstimate[i] + kData.offset_ptr()[i];
    loss += kData.weight_ptr()[i] * (kData.y_ptr()[i] * deltafunc_est -
                                     std::log(1.0 + std::exp(deltafunc_est)));
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

void CBernoulli::FitBestConstant(const CDataset& kData, const Bag& kBag,
                                 const double* kFuncEstimate,
                                 unsigned long num_terminalnodes,
                                 std::vector<double>& residuals,
                                 CCARTTree& tree) {
  unsigned long obs_num = 0;
  unsigned long node_num = 0;
  vector<double> numerator_vec(num_terminalnodes, 0.0);
  vector<double> denom_vec(num_terminalnodes, 0.0);

  for (obs_num = 0; obs_num < kData.get_trainsize(); obs_num++) {
    if (kBag.get_element(obs_num)) {
      numerator_vec[tree.get_node_assignments()[obs_num]] +=
          kData.weight_ptr()[obs_num] * residuals[obs_num];
      denom_vec[tree.get_node_assignments()[obs_num]] +=
          kData.weight_ptr()[obs_num] *
          (kData.y_ptr()[obs_num] - residuals[obs_num]) *
          (1 - kData.y_ptr()[obs_num] + residuals[obs_num]);
    }
  }

  for (node_num = 0; node_num < num_terminalnodes; node_num++) {
    if (tree.get_terminal_nodes()[node_num] != NULL) {
      if (denom_vec[node_num] == 0) {
        tree.get_terminal_nodes()[node_num]->set_prediction(0.0);
      } else {
        double temp = numerator_vec[node_num] / denom_vec[node_num];
        // avoid large changes in predictions on log odds scale
        if (std::abs(temp) > terminalnode_cap_level_) {
          if (!terminalnode_capped_) {
            // set fCappedPred=true so that warning only issued once
            terminalnode_capped_ = true;
            Rcpp::warning(
                "Some terminal node predictions were excessively large for "
                "Bernoulli and have been capped. Likely due to a "
                "feature that separates the 0/1 outcomes. Consider reducing "
                "shrinkage parameter.");
          }
          if (temp > terminalnode_cap_level_) {
            temp = terminalnode_cap_level_;
          } else if (temp < -terminalnode_cap_level_) {
            temp = -terminalnode_cap_level_;
          }
        }

        tree.get_terminal_nodes()[node_num]->set_prediction(temp);
      }
    }
  }
}

double CBernoulli::BagImprovement(const CDataset& kData, const Bag& kBag,
                                  const double* kFuncEstimate,
                                  const double kShrinkage,
                                  const std::vector<double>& kDeltaEstimate) {
  double returnvalue = 0.0;
  double deltafunc_est = 0.0;
  double weight = 0.0;
  unsigned long i = 0;

  for (i = 0; i < kData.get_trainsize(); i++) {
    if (!kBag.get_element(i)) {
      deltafunc_est = kFuncEstimate[i] + kData.offset_ptr()[i];

      if (kData.y_ptr()[i] == 1.0) {
        returnvalue += kData.weight_ptr()[i] * kShrinkage * kDeltaEstimate[i];
      }
      returnvalue += kData.weight_ptr()[i] *
                     (std::log(1.0 + std::exp(deltafunc_est)) -
                      std::log(1.0 + std::exp(deltafunc_est +
                                              kShrinkage * kDeltaEstimate[i])));
      weight += kData.weight_ptr()[i];
    }
  }

  return returnvalue / weight;
}
