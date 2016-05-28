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
// Description:
//
// Parameters:
//  radY           - R object (SEXP) containing the response of each data-point
//  radOffset      - R object (SEXP) containing the offset applied to each
//  prediction
//  radX           - R object (SEXP) containing the predictor values.
//  raiXOrder      - R object (SEXP) containing the order of predictor values to
//  be
//  used in GBM formula
//  radWeight      - R object (SEXP) containing weights to be used in fitting
//  process
//  racVarClasses  - R object (SEXP) containing the variable classes
//  ralMonotoneVar - R object (SEXP) containing +-1/9 indicating whether the
//  variables are:
//  monotone increasing (+1), decreasing (-1) or arbitrary (0) with the response
//  variable.
//	cTrain		   - int specifiy the number of data points in
//training set
//
//-----------------------------------
CDataset::CDataset(const DataDistParams& dataparams)
    : xmatrix_(dataparams.xvalues),
      response_(dataparams.response),
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

  // Set-up Bags
  databag_.assign(dataparams.num_trainrows, false);
  bagfraction_ = dataparams.bagfraction;
  totalinbag_ =
      (long)(dataparams.bagfraction * dataparams.num_trainobservations);

  // Ensure initialization makes sense
  if (totalinbag_ <= 0) {
    throw gbm_exception::InvalidArgument("you have an empty bag!");
  }
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
// Function: ~CDataset()
//
// Returns: none
//
// Description: default destructor for the dataset
//
// Parameters: none
//
//-----------------------------------
CDataset::~CDataset() {}

//-----------------------------------
// Function: random_order
//
// Returns: randomized order vector
//
// Parameters: none
//
//-----------------------------------
typedef std::vector<int> index_vector;
index_vector CDataset::RandomOrder() const {
  index_vector result(ncol());

  // fill the vector
  for (index_vector::size_type ind = 0; ind != result.size(); ++ind) {
    result[ind] = ind;
  }

  // and now shuffle
  std::random_shuffle(result.begin(), result.end(), gbm_functions::PtrShuffler);
  // and return
  return result;
}
