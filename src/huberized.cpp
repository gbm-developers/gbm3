//-----------------------------------
//
// File: huberized.cpp
//
// Description: huberized hinge loss.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "huberized.h"

//----------------------------------------
// Function Members - Private
//----------------------------------------
CHuberized::CHuberized() {}

//----------------------------------------
// Function Members - Public
//----------------------------------------

CDistribution* CHuberized::Create(DataDistParams& distparams) {
  return new CHuberized();
}

CHuberized::~CHuberized() {}

void CHuberized::ComputeWorkingResponse(const CDataset& kData, const Bag& kBag,
                                        const double* kFuncEstimate,
                                        std::vector<double>& residuals) {
  unsigned long i = 0;
  double delta_func_est = 0.0;

  for (i = 0; i < kData.get_trainsize(); i++) {
    delta_func_est = kFuncEstimate[i] + kData.offset_ptr()[i];
    if ((2 * kData.y_ptr()[i] - 1) * delta_func_est < -1) {
      residuals[i] = -4 * (2 * kData.y_ptr()[i] - 1);
    } else if (1 - (2 * kData.y_ptr()[i] - 1) * delta_func_est < 0) {
      residuals[i] = 0;
    } else {
      residuals[i] = -2 * (2 * kData.y_ptr()[i] - 1) *
                     (1 - (2 * kData.y_ptr()[i] - 1) * delta_func_est);
    }
  }
}

double CHuberized::InitF(const CDataset& kData) {
  unsigned long i = 0;
  double numerator = 0.0;
  double denominator = 0.0;

  for (i = 0; i < kData.get_trainsize(); i++) {
    if (kData.y_ptr()[i] == 1.0) {
      numerator += kData.weight_ptr()[i];
    } else {
      denominator += kData.weight_ptr()[i];
    }
  }

  return numerator / denominator;
}

double CHuberized::Deviance(const CDataset& kData, const Bag& kBag,
                            const double* kFuncEstimate) {
  unsigned long i = 0;
  double loss = 0.0;
  double delta_func_est = 0.0;
  double weights = 0.0;

  unsigned long num_rows_in_set = kData.get_size_of_set();

  for (i = 0; i < num_rows_in_set; i++) {
    delta_func_est = kData.offset_ptr()[i] + kFuncEstimate[i];
    if ((2 * kData.y_ptr()[i] - 1) * delta_func_est < -1) {
      loss += -kData.weight_ptr()[i] * 4 * (2 * kData.y_ptr()[i] - 1) *
              delta_func_est;
      weights += kData.weight_ptr()[i];
    } else if (1 - (2 * kData.y_ptr()[i] - 1) * delta_func_est < 0) {
      loss += 0;
      weights += kData.weight_ptr()[i];
    } else {
      loss += kData.weight_ptr()[i] *
              (1 - (2 * kData.y_ptr()[i] - 1) * delta_func_est) *
              (1 - (2 * kData.y_ptr()[i] - 1) * delta_func_est);
      weights += kData.weight_ptr()[i];
    }
  }  // close for(

  // TODO: Check if weights are all zero for validation set
  if ((weights == 0.0) && (loss == 0.0)) {
    return nan("");
  } else if (weights == 0.0) {
    return copysign(HUGE_VAL, loss);
  }

  return loss / weights;
}

void CHuberized::FitBestConstant(const CDataset& kData, const Bag& kBag,
                                 const double* kFuncEstimate,
                                 unsigned long num_terminalnodes,
                                 std::vector<double>& residuals,
                                 CCARTTree& tree) {
  double delta_func_est = 0.0;
  unsigned long obs_num = 0;
  unsigned long node_num = 0;

  vector<double> numerator_vec(num_terminalnodes, 0.0);
  vector<double> denominator_vec(num_terminalnodes, 0.0);

  for (obs_num = 0; obs_num < kData.get_trainsize(); obs_num++) {
    if (kBag.get_element(obs_num)) {
      delta_func_est = kFuncEstimate[obs_num] + kData.offset_ptr()[obs_num];
      if ((2 * kData.y_ptr()[obs_num] - 1) * kFuncEstimate[obs_num] < -1) {
        numerator_vec[tree.get_node_assignments()[obs_num]] +=
            kData.weight_ptr()[obs_num] * 4 * (2 * kData.y_ptr()[obs_num] - 1);
        denominator_vec[tree.get_node_assignments()[obs_num]] +=
            -kData.weight_ptr()[obs_num] * 4 *
            (2 * kData.y_ptr()[obs_num] - 1) * delta_func_est;
      } else if (1 - (2 * kData.y_ptr()[obs_num] - 1) * kFuncEstimate[obs_num] <
                 0) {
        numerator_vec[tree.get_node_assignments()[obs_num]] += 0;
        denominator_vec[tree.get_node_assignments()[obs_num]] += 0;
      } else {
        numerator_vec[tree.get_node_assignments()[obs_num]] +=
            kData.weight_ptr()[obs_num] * 2 * (2 * kData.y_ptr()[obs_num] - 1) *
            (1 - (2 * kData.y_ptr()[obs_num] - 1) * kFuncEstimate[obs_num]);
        denominator_vec[tree.get_node_assignments()[obs_num]] +=
            kData.weight_ptr()[obs_num] *
            (1 - (2 * kData.y_ptr()[obs_num] - 1) * kFuncEstimate[obs_num]) *
            (1 - (2 * kData.y_ptr()[obs_num] - 1) * kFuncEstimate[obs_num]);
      }
    }  // close if(afInBag[obs_num
  }

  for (node_num = 0; node_num < num_terminalnodes; node_num++) {
    if (tree.has_node(node_num)) {
      if (denominator_vec[node_num] == 0) {
        tree.get_terminal_nodes()[node_num]->set_prediction(0.0);
      } else {
        tree.get_terminal_nodes()[node_num]->set_prediction(
            numerator_vec[node_num] / denominator_vec[node_num]);
      }
    }
  }
}

double CHuberized::BagImprovement(const CDataset& kData, const Bag& kBag,
                                  const double* kFuncEstimate,
                                  const double kShrinkage,
                                  const std::vector<double>& kDeltaEstimate) {
  double returnvalue = 0.0;
  double delta_func_est = 0.0;
  double weight = 0.0;
  unsigned long i = 0;

  for (i = 0; i < kData.get_trainsize(); i++) {
    if (!kBag.get_element(i)) {
      delta_func_est = kFuncEstimate[i] + kData.offset_ptr()[i];

      if ((2 * kData.y_ptr()[i] - 1) * delta_func_est < -1) {
        returnvalue += kData.weight_ptr()[i] *
                       (-4 * (2 * kData.y_ptr()[i] - 1) * delta_func_est -
                        -4 * (2 * kData.y_ptr()[i] - 1) *
                            (delta_func_est + kShrinkage * kDeltaEstimate[i]));
        weight += kData.weight_ptr()[i];
      } else if (1 - (2 * kData.y_ptr()[i] - 1) * delta_func_est < 0) {
        returnvalue += 0;
        weight += kData.weight_ptr()[i];
      } else {
        returnvalue +=
            kData.weight_ptr()[i] *
            ((1 - (2 * kData.y_ptr()[i] - 1) * delta_func_est) *
                 (1 - (2 * kData.y_ptr()[i] - 1) * delta_func_est) -
             (1 -
              (2 * kData.y_ptr()[i] - 1) *
                  (delta_func_est + kShrinkage * kDeltaEstimate[i])) *
                 (1 -
                  (2 * kData.y_ptr()[i] - 1) *
                      (delta_func_est + kShrinkage * kDeltaEstimate[i])));
        // TODO: Does this require an weight+= ?
      }
    }
  }

  // TODO: Check if weights are all zero for validation set
  if ((weight == 0.0) && (returnvalue == 0.0)) {
    return nan("");
  } else if (weight == 0.0) {
    return copysign(HUGE_VAL, returnvalue);
  }

  return returnvalue / weight;
}
