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
CGaussian::CGaussian() {}

//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CGaussian::Create(DataDistParams& distparams) {
  return new CGaussian();
}

CGaussian::~CGaussian() {}

void CGaussian::ComputeWorkingResponse(const CDataset& kData,
                                       const double* kFuncEstimate,
                                       double* residuals) {
  unsigned long i = 0;

  if (!(kData.y_ptr() && kFuncEstimate && residuals && kData.weight_ptr())) {
    throw gbm_exception::InvalidArgument();
  }

  for (i = 0; i < kData.get_trainsize(); i++) {
    residuals[i] = kData.y_ptr()[i] - kData.offset_ptr()[i] - kFuncEstimate[i];
  }
}

double CGaussian::InitF(const CDataset& kData) {
  double sum = 0.0;
  double totalweight = 0.0;
  unsigned long i = 0;

  // compute the mean

  for (i = 0; i < kData.get_trainsize(); i++) {
    sum += kData.weight_ptr()[i] * (kData.y_ptr()[i] - kData.offset_ptr()[i]);
    totalweight += kData.weight_ptr()[i];
  }

  return sum / totalweight;
}

double CGaussian::Deviance(const CDataset& kData, const double* kFuncEstimate) {
  unsigned long i = 0;
  double loss = 0.0;
  double weight = 0.0;

  unsigned long num_rows_in_set = kData.get_size_of_set();
  for (i = 0; i < num_rows_in_set; i++) {
    loss += kData.weight_ptr()[i] *
            (kData.y_ptr()[i] - kData.offset_ptr()[i] - kFuncEstimate[i]) *
            (kData.y_ptr()[i] - kData.offset_ptr()[i] - kFuncEstimate[i]);
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

void CGaussian::FitBestConstant(const CDataset& kData,
                                const double* kFuncEstimate,
                                unsigned long num_terminalnodes,
                                double* residuals, CCARTTree& tree) {
  // the tree aready stores the mean prediction
  // no refitting necessary
}

double CGaussian::BagImprovement(const CDataset& kData,
                                 const double* kFuncEstimate,
                                 const double kShrinkage,
                                 const double* kDeltaEstimate) {
  double returnvalue = 0.0;
  double deltafunc_est = 0.0;
  double weight = 0.0;
  unsigned long i = 0;

  for (i = 0; i < kData.get_trainsize(); i++) {
    if (!kData.get_bag_element(i)) {
      deltafunc_est = kFuncEstimate[i] + kData.offset_ptr()[i];

      returnvalue += kData.weight_ptr()[i] * kShrinkage * kDeltaEstimate[i] *
                     (2.0 * (kData.y_ptr()[i] - deltafunc_est) -
                      kShrinkage * kDeltaEstimate[i]);
      weight += kData.weight_ptr()[i];
    }
  }

  return returnvalue / weight;
}
