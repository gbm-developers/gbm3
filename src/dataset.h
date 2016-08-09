//------------------------------------------------------------------------------
//
//  File:       dataset.h
//
//  Description:   Header file for dataset public methods.
//
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//				31/03/2016  James Hickey: RAII and Pimpled
//				25/05/2016  James Hickey: Strip out pimpl
//------------------------------------------------------------------------------

#ifndef DATASET_H
#define DATASET_H

//------------------------------
// Includes
//------------------------------
#include "datadistparams.h"
#include "gbm_exception.h"
#include "gbm_functions.h"
#include <algorithm>
#include <memory>
#include <vector>
#include <Rcpp.h>

//------------------------------
// Type Defs
//------------------------------
typedef std::vector<int> index_vector;

//------------------------------
// Class definition
//------------------------------
class CDataset {
 public:
  //----------------------
  // Public Constructors
  //----------------------
  CDataset(const DataDistParams& dataParams);

  //---------------------
  // Public Functions
  //---------------------
  unsigned int nrow() const { return xmatrix_.nrow(); };
  unsigned int ncol() const { return xmatrix_.ncol(); };

  double* y_ptr(long colIndex = 0) {
    return yptrs_[colIndex];
  };  // get iterator to class labels

  const double* y_ptr(long colIndex = 0) const {
    return yptrs_[colIndex];
  };  // const overloaded version

  int* yint_ptr(long colIndex = 0) {
    return yintptrs_[colIndex];
  };  // get iterator to class labels

  const int* yint_ptr(long colIndex = 0) const {
    return yintptrs_[colIndex];
  };  // const overloaded version

  const double* offset_ptr() const { return offset_ptr_; };

  const double* weight_ptr() const {
    return weights_ptr_;
  };  // const overloaded version

  int varclass(int ind) const { return num_variable_classes_[ind]; };
  int monotone(int ind) const { return variable_monotonicity_[ind]; };

  const int* order_ptr() const { return order_xvals_.begin(); };

  double x_value(const int row, const int col) const {
    return xmatrix_(row, col);
  };  // retrieve predictor value

  unsigned long get_trainsize() const {
    return num_traindata_;
  };  // get size of training set
  unsigned long get_num_features() const {
    return num_features_;
  };  // get the number of features in data

  void shift_to_validation() {
    if (point_at_trainingset_) {
      for (unsigned int i = 0; i < yptrs_.size(); i++) {
        yptrs_[i] = shift_ptr_to_validation(yptrs_[i]);
      }
      for (unsigned int i = 0; i < yintptrs_.size(); i++) {
        yintptrs_[i] = shift_ptr_to_validation(yintptrs_[i]);
      }
      offset_ptr_ = shift_ptr_to_validation(offset_ptr_);
      weights_ptr_ = shift_ptr_to_validation(weights_ptr_);
      point_at_trainingset_ = false;
    } else {
      throw gbm_exception::InvalidArgument(
          "Data is already the validation set.");
    }
  };
  void shift_to_train() {
    if (!(point_at_trainingset_)) {
      for (unsigned int i = 0; i < yptrs_.size(); i++) {
        yptrs_[i] = shift_ptr_to_train(yptrs_[i]);
      }
      for (unsigned int i = 0; i < yintptrs_.size(); i++) {
        yintptrs_[i] = shift_ptr_to_train(yintptrs_[i]);
      }
      offset_ptr_ = shift_ptr_to_train(offset_ptr_);
      weights_ptr_ = shift_ptr_to_train(weights_ptr_);
      point_at_trainingset_ = true;
    } else {
      throw gbm_exception::InvalidArgument("Data is already the training set.");
    }
  };

  index_vector RandomOrder() const;  // randomize order of predictor varaiables
  unsigned long get_validsize() const { return num_validationdata_; };
  unsigned long get_size_of_set() const {
    if (point_at_trainingset_) return get_trainsize();
    return get_validsize();
  }
  unsigned long get_num_observations_in_training() const {
    return num_trainobservations_;
  }
  int get_row_observation_id(int row_number) const {
    return observation_ids_(row_number);
  }

 private:
  //-------------------
  // Private Variables
  //-------------------
  void set_up_yptrs() {
    for (long i = 0; i < response_.ncol(); i++) {
      yptrs_.push_back(response_(Rcpp::_, i).begin());
    }
    for (long i = 0; i < intResponse_.ncol(); i++) {
      yintptrs_.push_back(intResponse_(Rcpp::_, i).begin());
    }
  }

  template <typename T>
  T* shift_ptr_to_validation(T* x) const {
    if (x) {
      return x + num_traindata_;
    } else {
      return x;
    }
  }

  template <typename T>
  T* shift_ptr_to_train(T* x) const {
    if (x) {
      return x - num_traindata_;
    } else {
      return x;
    }
  }
  //-------------------
  // Private Variables
  //-------------------

  // Numeric vectors storing data
  Rcpp::NumericMatrix xmatrix_, response_;
  Rcpp::IntegerMatrix intResponse_;
  Rcpp::NumericVector response_offset_, data_weights_;
  Rcpp::IntegerVector num_variable_classes_, variable_monotonicity_,
      order_xvals_, observation_ids_;

  // Ptrs to numeric vectors - these must be mutable
  std::vector<double*> yptrs_;
  std::vector<int*> yintptrs_;
  double* offset_ptr_;
  double* weights_ptr_;

  // Properties of the data
  unsigned long num_traindata_;
  unsigned long num_trainobservations_;
  unsigned long num_validationdata_;
  unsigned long num_features_;
  bool point_at_trainingset_;
};
#endif  // DATASET_H
