//-----------------------------------
//
// File: tdist.cpp
//
// Description: t distribution implementation for GBM.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "locationm.h"
#include "tdist.h"
#include <vector>

//----------------------------------------
// Function Members - Private
//----------------------------------------
CTDist::CTDist(double nu) : mplocm_("tdist", nu) { m_nu_ = nu; }

//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CTDist::Create(DataDistParams& distparams) {
  // Check that misc exists
  double nu = Rcpp::as<double>(distparams.misc[0]);
  if (!gbm_functions::has_value(nu)) {
    throw gbm_exception::Failure("T Dist requires misc to initialization.");
  }
  return new CTDist(nu);
}

CTDist::~CTDist() {}

void CTDist::ComputeWorkingResponse(const CDataset& kData,
                                    const double* kFuncEstimate,
                                    std::vector<double>& residuals) {
  unsigned long i = 0;
  double du = 0.0;

  for (i = 0; i < kData.get_trainsize(); i++) {
    du = kData.y_ptr()[i] - kData.offset_ptr()[i] - kFuncEstimate[i];
    residuals[i] = (2 * du) / (m_nu_ + (du * du));
  }
}

double CTDist::InitF(const CDataset& kData) {
  // Get objects to pass into the LocM function
  std::vector<double> arr(kData.get_trainsize());

  for (unsigned long ii = 0; ii < kData.get_trainsize(); ii++) {
    double offset = kData.offset_ptr()[ii];
    arr[ii] = kData.y_ptr()[ii] - offset;
  }

  return mplocm_.LocationM(kData.get_trainsize(), &arr[0], kData.weight_ptr(),
                           0.5);
}

double CTDist::Deviance(const CDataset& kData, const double* kFuncEstimate) {
  unsigned long i = 0;
  double loss = 0.0;
  double weight = 0.0;
  double du = 0.0;

  // Switch to validation set if necessary
  unsigned long num_rows_in_set = kData.get_size_of_set();

  for (i = 0; i < num_rows_in_set; i++) {
    du = kData.y_ptr()[i] - kData.offset_ptr()[i] - kFuncEstimate[i];
    loss += kData.weight_ptr()[i] * std::log(m_nu_ + (du * du));
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

void CTDist::FitBestConstant(const CDataset& kData, const double* kFuncEstimate,
                             unsigned long num_terminalnodes, std::vector<double>& residuals,
                             CCARTTree& tree) {
  // Local variables
  unsigned long node_num = 0;
  unsigned long obs_num = 0;

  std::vector<double> arr_vec, weight_vec;
  // Call LocM for the array of values on each node
  for (node_num = 0; node_num < num_terminalnodes; node_num++) {
    if (tree.get_terminal_nodes()[node_num]->get_numobs() >=
        tree.min_num_obs_required()) {
      arr_vec.clear();
      weight_vec.clear();

      for (obs_num = 0; obs_num < kData.get_trainsize(); obs_num++) {
        if (kData.get_bag_element(obs_num) &&
            (tree.get_node_assignments()[obs_num] == node_num)) {
          const double dOffset = kData.offset_ptr()[obs_num];
          arr_vec.push_back(kData.y_ptr()[obs_num] - dOffset -
                            kFuncEstimate[obs_num]);
          weight_vec.push_back(kData.weight_ptr()[obs_num]);
        }
      }

      tree.get_terminal_nodes()[node_num]->set_prediction(
          mplocm_.LocationM(arr_vec.size(), &arr_vec[0], &weight_vec[0], 0.5));
    }
  }
}

double CTDist::BagImprovement(const CDataset& kData,
                              const double* kFuncEstimate,
                              const double kShrinkage,
                              const double* kDeltaEstimate) {
  double returnvalue = 0.0;
  unsigned long i = 0;
  double weight = 0.0;

  for (i = 0; i < kData.get_trainsize(); i++) {
    if (!kData.get_bag_element(i)) {
      const double dF = kFuncEstimate[i] + kData.offset_ptr()[i];
      const double dU = (kData.y_ptr()[i] - dF);
      const double dV =
          (kData.y_ptr()[i] - dF - kShrinkage * kDeltaEstimate[i]);

      returnvalue += kData.weight_ptr()[i] *
                     (std::log(m_nu_ + (dU * dU)) - log(m_nu_ + (dV * dV)));
      weight += kData.weight_ptr()[i];
    }
  }

  return returnvalue / weight;
}
