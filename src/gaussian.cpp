//-----------------------------------
//
// File: gaussian.cpp
//
// Description: gaussian distribution implementation for GBM.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "gaussian.h"

//----------------------------------------
// Function Members - Private
//----------------------------------------
CGaussian::CGaussian(const parallel_details& parallel)
    : CDistribution(parallel) {}

//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CGaussian::Create(DataDistParams& distparams) {
  return new CGaussian(distparams.parallel);
}

CGaussian::~CGaussian() {}

void CGaussian::ComputeWorkingResponse(const CDataset& kData, const Bag& kBag,
                                       const double* kFuncEstimate,
                                       std::vector<double>& residuals) {
  if (!(kData.y_ptr() && kFuncEstimate &&
        kData.weight_ptr())) {
    throw gbm_exception::InvalidArgument();
  }

#pragma omp parallel for schedule(static, get_array_chunk_size()) \
  num_threads(get_num_threads())
  for (unsigned long i = 0; i < kData.get_trainsize(); i++) {
    residuals[i] = kData.y_ptr()[i] - kData.offset_ptr()[i] - kFuncEstimate[i];
  }
}

double CGaussian::InitF(const CDataset& kData) {
  double sum = 0.0;
  double totalweight = 0.0;

// compute the mean

#pragma omp parallel for schedule(static, get_array_chunk_size()) \
    reduction(+ : sum, totalweight) num_threads(get_num_threads())
  for (unsigned long i = 0; i < kData.get_trainsize(); i++) {
    sum += kData.weight_ptr()[i] * (kData.y_ptr()[i] - kData.offset_ptr()[i]);
    totalweight += kData.weight_ptr()[i];
  }

  return sum / totalweight;
}

double CGaussian::Deviance(const CDataset& kData, const Bag& kBag,
                           const double* kFuncEstimate) {
  double loss = 0.0;
  double weight = 0.0;

  unsigned long num_rows_in_set = kData.get_size_of_set();
#pragma omp parallel for schedule(static, get_array_chunk_size()) \
    reduction(+ : loss, weight) num_threads(get_num_threads())
  for (unsigned long i = 0; i < num_rows_in_set; i++) {
    const double tmp =
        (kData.y_ptr()[i] - kData.offset_ptr()[i] - kFuncEstimate[i]);
    loss += kData.weight_ptr()[i] * tmp * tmp;
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

void CGaussian::FitBestConstant(const CDataset& kData, const Bag& kBag,
                                const double* kFuncEstimate,
                                unsigned long num_terminalnodes,
                                std::vector<double>& residuals,
                                CCARTTree& tree) {
  // the tree aready stores the mean prediction
  // no refitting necessary
}

double CGaussian::BagImprovement(const CDataset& kData, const Bag& kBag,
                                 const double* kFuncEstimate,
                                 const double kShrinkage,
                                 const std::vector<double>& kDeltaEstimate) {
  double returnvalue = 0.0;
  double weight = 0.0;

#pragma omp parallel for schedule(static, get_array_chunk_size()) \
    reduction(+ : returnvalue, weight) num_threads(get_num_threads())
  for (unsigned long i = 0; i < kData.get_trainsize(); i++) {
    if (!kBag.get_element(i)) {
      const double deltafunc_est = kFuncEstimate[i] + kData.offset_ptr()[i];

      returnvalue += kData.weight_ptr()[i] * kShrinkage * kDeltaEstimate[i] *
                     (2.0 * (kData.y_ptr()[i] - deltafunc_est) -
                      kShrinkage * kDeltaEstimate[i]);
      weight += kData.weight_ptr()[i];
    }
  }

  return returnvalue / weight;
}
