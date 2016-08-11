//-----------------------------------
//
// File: quantile.cpp
//
// Description: quantile distribution for GBM.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "quantile.h"

//----------------------------------------
// Function Members - Private
//----------------------------------------
CQuantile::CQuantile(double alpha, const parallel_details& parallel)
    : CDistribution(parallel), mplocm_("Other") {
  alpha_ = alpha;
}

//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CQuantile::Create(DataDistParams& distparams) {
  // Check that misc exists
  // This isn't very nice but helps with garbage collection
  // - only one Rcpp member pointing at a SEXP
  Rcpp::List misc_list = Rcpp::as<Rcpp::List>(distparams.misc);
  double alpha = Rcpp::as<double>(misc_list[0]);
  if (!gbm_functions::has_value(alpha)) {
    throw gbm_exception::Failure(
        "Quantile dist requires misc to initialization.");
  }
  return new CQuantile(alpha, distparams.parallel);
}

CQuantile::~CQuantile() {}

void CQuantile::ComputeWorkingResponse(const CDataset& kData, const Bag& kBag,
                                       const double* kFuncEstimate,
                                       std::vector<double>& residuals) {
#pragma omp parallel for schedule(static, get_array_chunk_size()) \
  num_threads(get_num_threads())
  for (unsigned long i = 0; i < kData.get_trainsize(); i++) {
    residuals[i] = (kData.y_ptr()[i] > kFuncEstimate[i] + kData.offset_ptr()[i])
                       ? alpha_
                       : -(1.0 - alpha_);
  }
}

double CQuantile::InitF(const CDataset& kData) {
  vecd_.resize(kData.get_trainsize());
#pragma omp parallel for schedule(static, get_array_chunk_size()) \
  num_threads(get_num_threads())
  for (unsigned long i = 0; i < kData.get_trainsize(); i++) {
    vecd_[i] = kData.y_ptr()[i] - kData.offset_ptr()[i];
  }

  return mplocm_.WeightedQuantile(kData.get_trainsize(), &vecd_[0],
                                  kData.weight_ptr(), alpha_);
}

double CQuantile::Deviance(const CDataset& kData, const Bag& kBag,
                           const double* kFuncEstimate) {
  double loss = 0.0;
  double weight = 0.0;

  // Switch to validation set if necessary
  unsigned long num_rows_in_set = kData.get_size_of_set();

#pragma omp parallel for schedule(static, get_array_chunk_size()) \
    reduction(+ : loss, weight) num_threads(get_num_threads())
  for (unsigned long i = 0; i < num_rows_in_set; i++) {
    if (kData.y_ptr()[i] > kFuncEstimate[i] + kData.offset_ptr()[i]) {
      loss += kData.weight_ptr()[i] * alpha_ *
              (kData.y_ptr()[i] - kFuncEstimate[i] - kData.offset_ptr()[i]);
    } else {
      loss += kData.weight_ptr()[i] * (1.0 - alpha_) *
              (kFuncEstimate[i] + kData.offset_ptr()[i] - kData.y_ptr()[i]);
    }
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

void CQuantile::FitBestConstant(const CDataset& kData, const Bag& kBag,
                                const double* kFuncEstimate,
                                unsigned long num_terminalnodes,
                                std::vector<double>& residuals,
                                CCARTTree& tree) {
  unsigned long node_num = 0;
  unsigned long obs_num = 0;
  unsigned long vec_num = 0;

  vecd_.resize(
      kData.get_trainsize());  // should already be this size from InitF
  std::vector<double> weight_vec(kData.get_trainsize());

  for (node_num = 0; node_num < num_terminalnodes; node_num++) {
    if (tree.get_terminal_nodes()[node_num]->get_numobs() >=
        tree.min_num_obs_required()) {
      vec_num = 0;
      for (obs_num = 0; obs_num < kData.get_trainsize(); obs_num++) {
        if (kBag.get_element(obs_num) &&
            (tree.get_node_assignments()[obs_num] == node_num)) {
          const double offset = kData.offset_ptr()[obs_num];

          vecd_[vec_num] =
              kData.y_ptr()[obs_num] - offset - kFuncEstimate[obs_num];
          weight_vec[vec_num] = kData.weight_ptr()[obs_num];
          vec_num++;
        }
      }

      tree.get_terminal_nodes()[node_num]->set_prediction(
          mplocm_.WeightedQuantile(vec_num, &vecd_[0], &weight_vec[0], alpha_));
    }
  }
}

double CQuantile::BagImprovement(const CDataset& kData, const Bag& kBag,
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

      if (kData.y_ptr()[i] > delta_func_est) {
        returnvalue += kData.weight_ptr()[i] * alpha_ *
                       (kData.y_ptr()[i] - delta_func_est);
      } else {
        returnvalue += kData.weight_ptr()[i] * (1 - alpha_) *
                       (delta_func_est - kData.y_ptr()[i]);
      }

      if (kData.y_ptr()[i] > delta_func_est + kShrinkage * kDeltaEstimate[i]) {
        returnvalue -= kData.weight_ptr()[i] * alpha_ *
                       (kData.y_ptr()[i] - delta_func_est -
                        kShrinkage * kDeltaEstimate[i]);
      } else {
        returnvalue -= kData.weight_ptr()[i] * (1 - alpha_) *
                       (delta_func_est + kShrinkage * kDeltaEstimate[i] -
                        kData.y_ptr()[i]);
      }
      weight += kData.weight_ptr()[i];
    }
  }
  return returnvalue / weight;
}
