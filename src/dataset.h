//------------------------------------------------------------------------------
//
//  File:       dataset.h
//
//  Description:   Header file for dataset public methods. The dataset
//    accessed via a pointer to its implementation.
//
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//				31/03/2016  James Hickey: RAII and Pimpled
//------------------------------------------------------------------------------

#ifndef DATASET_H
#define DATASET_H

//------------------------------
// Includes
//------------------------------
#include "buildinfo.h"
#include "config_structs.h"
#include "gbmexcept.h"
#include "gbm_functions.h"
#include <algorithm>
#include <memory>
#include <vector>
#include <Rcpp.h>

typedef std::vector<int> bag;

//-----------------------------------
// Class Definition - Private Variable
//-----------------------------------
class CDImpl
{
public:
  //----------------------
  // Public Constructors
  //----------------------
  CDImpl(SEXP response, SEXP offset_vec, SEXP covariates, SEXP covar_order,
	 SEXP obs_weight, SEXP var_classes, SEXP monotonicity_vec,
	 int num_trainrows, const int num_features , double fraction_inbag,
	 int num_trainobs, SEXP observationids) :
    xmatrix_(covariates), response_(response), response_offset_(offset_vec), data_weights_(obs_weight),
    num_variable_classes_(var_classes), variable_monotonicity_(monotonicity_vec),
    order_xvals_(covar_order), patient_ids_(observationids)
  {
    
    // If you've no offset set to 0
    if(!GBM_FUNC::has_value(response_offset_))
    {
      Rcpp::NumericVector new_offset(xmatrix_.nrow());
      std::swap(response_offset_, new_offset);
    }
    
    // Set-up pointers
    set_up_yptrs();
    weights_ptr_ = data_weights_.begin();
    offset_ptr_ = response_offset_.begin();
    
    // Set-up data properties
    num_traindata_ = num_trainrows;
    num_trainobservations_ = num_trainobs;
    num_validationdata_ = xmatrix_.nrow() - num_trainrows;
    num_features_  = num_features;
    point_at_trainingset_ = true;
    
    // Set-up Bags
    databag_.assign(num_trainrows, false);
    bagfraction_ = fraction_inbag;
    totalinbag_ = (long) (fraction_inbag * num_trainobs);

    // Ensure initialization makes sense
    if (totalinbag_ <= 0)
      {
	throw GBM::InvalidArgument("you have an empty bag!");
      }
    if (num_trainrows <= 0)
      {
	throw GBM::InvalidArgument("you've <= 0 training instances");
      }
    
  };
  
  //---------------------
  // Public destructor
  //---------------------
  ~CDImpl(){};
  
  //---------------------
  // Public Functions
  //---------------------
  //-----------------------------------
  // Function: shift_ptr_to_validation
  //
  // Returns:  shifts the ptr to the validation set.
  //
  // Parameters: none
  //
  //-----------------------------------
  template<typename T>
    T* shift_ptr_to_validation(T* x) const
    {
      if(x)
      {
    	  return x + num_traindata_;
      }
      else
      {
    	  return x;
      }
    }
  
  //-----------------------------------
  // Function: shift_ptr_to_train
  //
  // Returns:  shifts the ptr to the training set.
  //
  // Parameters: none
  //
  //-----------------------------------
  template<typename T>
	T* shift_ptr_to_train(T* x) const
	{
	  if(x)
	  {
		  return x - num_traindata_;
	  }
	  else
	  {
		  return x;
	  }
	}
  
  //-----------------------------------
  // Function: SetUpYPtrs
  //
  // Returns:  sets up the ptrs to each column of response mat.
  //
  // Parameters: none
  //
  //-----------------------------------
  void set_up_yptrs()
  {
    for(long i = 0; i < response_.ncol(); i++)
      {
    	yptrs_.push_back(response_(Rcpp::_, i).begin());
      }
  }
  
  void shift_to_train() const
  {
    if(!(point_at_trainingset_))
    {
      for(unsigned int i = 0; i < yptrs_.size(); i++)
      {
    	  yptrs_[i] = shift_ptr_to_train(yptrs_[i]);
      }
      offset_ptr_ = shift_ptr_to_train(offset_ptr_);
      weights_ptr_ = shift_ptr_to_train(weights_ptr_);
      point_at_trainingset_ = true;
    }
    else
    {
      throw GBM::InvalidArgument("Data is already the training set.");
    }
  }

  //---shift_to_validation() const
  void shift_to_validation() const
  {
    if (point_at_trainingset_)
    {
      for(unsigned int i = 0; i < yptrs_.size(); i++)
      {
    	  yptrs_[i] = shift_ptr_to_validation(yptrs_[i]);
      }
      offset_ptr_ = shift_ptr_to_validation(offset_ptr_);
      weights_ptr_ = shift_ptr_to_validation(weights_ptr_);
      point_at_trainingset_ = false;
    }
    else
    {
    	throw GBM::InvalidArgument("Data is already the validation set.");
    }
  }

  // Public Variables
  //-------------------
  // Numeric vectors storing data
  Rcpp::NumericMatrix xmatrix_, response_;
  Rcpp::NumericVector response_offset_, data_weights_;
  Rcpp::IntegerVector num_variable_classes_, variable_monotonicity_, order_xvals_, patient_ids_;
  
  // Ptrs to numeric vectors - these must be mutable
  mutable std::vector<double*> yptrs_;
  mutable double* offset_ptr_;
  mutable double* weights_ptr_;
  
  // Properties of the data
  unsigned long num_traindata_;
  long num_trainobservations_;
  unsigned long num_validationdata_;
  long num_features_;
  mutable bool point_at_trainingset_;
  
  // Bagged  data
  bag databag_;
  double bagfraction_;
  unsigned long totalinbag_;
  
};

//------------------------------
// Class definition
//------------------------------
class CDataset
{
public:

  //----------------------
  // Public Constructors
  //----------------------
  CDataset(const DataDistParams& dataParams);
  
  //---------------------
  // Public destructor
  //---------------------
  virtual ~CDataset();
  
  //---------------------
  // Public Functions
  //---------------------
  unsigned int nrow() const
  {
    return dataImpl.xmatrix_.nrow();
  };
  unsigned int ncol() const
  {
    return dataImpl.xmatrix_.ncol();
  };
  
  double* y_ptr(long colIndex=0)
  {
    return dataImpl.yptrs_[colIndex];
  }; //get iterator to class labels

  const double* y_ptr(long colIndex=0) const
  {
    return dataImpl.yptrs_[colIndex];
  }; //const overloaded version
  
  const double* offset_ptr() const
  {
    return dataImpl.offset_ptr_;
  };
  
  const double* weight_ptr() const
  {
    return dataImpl.weights_ptr_;
  }; //const overloaded version
  
  int varclass(int ind) const
  {
    return dataImpl.num_variable_classes_[ind];
  };
  int monotone(int ind) const
  {
    return dataImpl.variable_monotonicity_[ind];
  };
  
  const int* order_ptr() const
  {
    return dataImpl.order_xvals_.begin();
  };

  double x_value(const int row, const int col) const
  {
    return dataImpl.xmatrix_(row, col);
  }; // retrieve predictor value
  
  unsigned long get_trainsize() const { return dataImpl.num_traindata_; }; // get size of training set
  long get_num_features() const { return dataImpl.num_features_; }; // get the number of features in data
  
  void shift_to_validation() const {  dataImpl.shift_to_validation(); }; // shift all of the ptrs to validation set
  void shift_to_train() const {  dataImpl.shift_to_train(); }; // shift all of the ptrs to training set
  
  typedef std::vector<int> index_vector;
  index_vector RandomOrder() const;//randomize order of predictor varaiables
  
  double get_bagfraction() const { return dataImpl.bagfraction_; };
  
  unsigned long get_validsize() const { return dataImpl.num_validationdata_; };
  unsigned long get_total_in_bag() const { return dataImpl.totalinbag_;};
  
  unsigned long get_num_patients_in_training() const
  {
	  return dataImpl.num_trainobservations_;
  }

  int get_row_patient_id(int rowNumber) const
  {
	  return dataImpl.patient_ids_(rowNumber);
  }

  bool get_bag_element(long index) const
  {
	return dataImpl.databag_[index];
  }

  void set_bag_element(long index) { dataImpl.databag_[index] = 1; };

  void clear_bag()
  {
    dataImpl.databag_.assign(get_trainsize(), 0);
  };
  
private:
  //-------------------
  // Private Variables
  //-------------------
  CDImpl dataImpl;
};
#endif // DATASET_H


