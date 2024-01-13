//-----------------------------------
//
// File: dataset.cpp
//
// Description: class that implements the public methods and pimpl
//    of the dataset.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "dataset.h"
#include <Rcpp.h>

//----------------------------------------
// Function Members - Public
//----------------------------------------
//-----------------------------------
// Function: CDataset
//
// Returns: none
//
// Description: Constructs a d
//
// Parameters:
//  dataparams - a struct containing all of the data parameters
//				required to make a dataset object.
//-----------------------------------
CDataset::CDataset(const DataDistParams& dataparams)
    : xmatrix_(dataparams.xvalues),
      response_(dataparams.response),
      intResponse_(dataparams.intResponse),
      response_offset_(dataparams.offset),
      data_weights_(dataparams.variable_weight),
      num_variable_classes_(dataparams.variable_num_classes),
      variable_monotonicity_(dataparams.variable_monotonicity),
      order_xvals_(dataparams.xorder),
      observation_ids_(dataparams.observationids) {
  // If you've no offset set to 0
  if (!gbm_functions::has_value(response_offset_)) {
    Rcpp::NumericVector new_offset(xmatrix_.nrow());
    std::swap(response_offset_, new_offset);
  }

  // Set-up pointers
  set_up_yptrs();
  weights_ptr_ = data_weights_.begin();
  offset_ptr_ = response_offset_.begin();

  // Set-up data properties
  num_traindata_ = dataparams.num_trainrows;
  num_trainobservations_ = dataparams.num_trainobservations;
  num_validationdata_ = xmatrix_.nrow() - dataparams.num_trainrows;
  num_features_ = dataparams.num_features;
  point_at_trainingset_ = true;

  // Ensure initialization makes sense
  if (num_traindata_ <= 0) {
    throw gbm_exception::InvalidArgument("you've <= 0 training instances");
  }
  // Check for errors on initialization
  if (xmatrix_.ncol() != variable_monotonicity_.size()) {
    throw gbm_exception::InvalidArgument(
        "shape mismatch (monotone does not match data)");
  }

  if (xmatrix_.ncol() != num_variable_classes_.size()) {
    throw gbm_exception::InvalidArgument(
        "shape mismatch (var classes does not match data)");
  }

  if (xmatrix_.nrow() < int(dataparams.num_trainrows)) {
    throw gbm_exception::InvalidArgument(
        "your training instances don't make sense");
  }
}

//-----------------------------------
// Function: random_order
//
// Returns: randomized order vector
//
// Parameters: none
//
//-----------------------------------
index_vector CDataset::RandomOrder() const {
  index_vector result(ncol());

  // fill the vector
  for (index_vector::size_type ind = 0; ind != result.size(); ++ind) {
    result[ind] = ind;
  }

  // and now shuffle
  // std::random_shuffle(result.begin(), result.end(), gbm_functions::PtrShuffler);
  std::shuffle(result.begin(), result.end(), std::default_random_engine());
  // and return
  return result;
}
