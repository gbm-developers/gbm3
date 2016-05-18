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

#ifndef __dataset_h__
#define __dataset_h__

//------------------------------
// Includes
//------------------------------
#include "buildinfo.h"
#include "configStructs.h"
#include "gbmexcept.h"
#include "gbmFunc.h"
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
  CDImpl(SEXP radY, SEXP radOffset, SEXP radX, SEXP raiXOrder,
	 SEXP radWeight, SEXP racVarClasses, SEXP ralMonotoneVar,
	 const int cTrain, const int cFeatures, const double fractionInBag) :
    adY(radY), adOffset(radOffset), adWeight(radWeight), adX(radX),
    acVarClasses(racVarClasses), alMonotoneVar(ralMonotoneVar),
    aiXOrder(raiXOrder), numOfTrainData(cTrain), numOfFeatures(cFeatures)
  {
    
    // If you've no offset set to 0
    if(!GBM_FUNC::has_value(adOffset))
    {
      Rcpp::NumericVector new_offset(adX.nrow());
      std::swap(adOffset, new_offset);
    }
    
    // Set variables
    bagFraction = fractionInBag;
    totalInBag = (long) (fractionInBag * cTrain);
    cValid = adX.nrow() - cTrain;
    pointAtTrainSet = true;
    
    // Set up pointers
    adYPtr = adY(Rcpp::_, 0).begin();
    adWeightPtr = adWeight.begin();
    adOffsetPtr = adOffset.begin();
    afInBag.assign(cTrain, false);
    SetUpYPtrs();
    
    // Ensure initialization makes sense
    if (totalInBag <= 0)
      {
	throw GBM::invalid_argument("you have an empty bag!");
      }
    if (cTrain <= 0)
      {
	throw GBM::invalid_argument("you've <= 0 training instances");
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
    	  return x + numOfTrainData;
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
		  return x - numOfTrainData;
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
  void SetUpYPtrs()
  {
    for(long i = 0; i < adY.ncol(); i++)
      {
    	yptrs.push_back(adY(Rcpp::_, i).begin());
      }
  }
  
  void shift_to_train() const
  {
    if(!(pointAtTrainSet))
    {
      for(int i = 0; i < yptrs.size(); i++)
      {
    	  yptrs[i] = shift_ptr_to_train(yptrs[i]);
      }
      adOffsetPtr = shift_ptr_to_train(adOffsetPtr);
      adWeightPtr = shift_ptr_to_train(adWeightPtr);
      pointAtTrainSet = true;
    }
    else
    {
      throw GBM::invalid_argument("Data is already the training set.");
    }
  }

  //---shift_to_validation() const
  void shift_to_validation() const
  {
    if (pointAtTrainSet)
    {
      for(int i = 0; i < yptrs.size(); i++)
      {
    	  yptrs[i] = shift_ptr_to_validation(yptrs[i]);
      }
      adOffsetPtr = shift_ptr_to_validation(adOffsetPtr);
      adWeightPtr = shift_ptr_to_validation(adWeightPtr);
      pointAtTrainSet = false;
    }
    else
    {
    	throw GBM::invalid_argument("Data is already the validation set.");
    }
  }

  // Public Variables
  //-------------------
  // Numeric vectors storing data
  Rcpp::NumericVector adOffset, adWeight;
  Rcpp::NumericMatrix adX, adY;
  Rcpp::IntegerVector acVarClasses, alMonotoneVar, aiXOrder;
  
  // Ptrs to numeric vectors - these must be mutable
  mutable std::vector<double*> yptrs;
  mutable double* adYPtr;
  mutable double* adOffsetPtr;
  mutable double* adWeightPtr;
  
  // Properties of the data
  long numOfTrainData;
  unsigned long cValid;
  long numOfFeatures;
  mutable bool pointAtTrainSet;
  
  // Bagged  data
  bag afInBag;
  double bagFraction;
  long totalInBag;
  
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
  int nrow() const
  {
    return dataImpl.adX.nrow();
  };
  int ncol() const
  {
    return dataImpl.adX.ncol();
  };
  
  double* y_ptr(long colIndex=0)
  {
    return dataImpl.yptrs[colIndex];
  }; //get iterator to class labels

  const double* y_ptr(long colIndex=0) const
  {
    return dataImpl.yptrs[colIndex];
  }; //const overloaded version
  
  const double* offset_ptr() const
  {
    return dataImpl.adOffsetPtr;
  };
  
  const double* weight_ptr() const
  {
    return dataImpl.adWeightPtr;
  }; //const overloaded version
  
  int varclass(int ind) const
  {
    return dataImpl.acVarClasses[ind];
  };
  int monotone(int ind) const
  {
    return dataImpl.alMonotoneVar[ind];
  };
  
  const int* order_ptr() const
  {
    return dataImpl.aiXOrder.begin();
  };

  double x_value(const int row, const int col) const
  {
    return dataImpl.adX(row, col);
  }; // retrieve predictor value
  
  long get_trainSize() const { return dataImpl.numOfTrainData; }; // get size of training set
  long get_numFeatures() const; // get the number of features in data
  
  void shift_to_validation() const; // shift all of the ptrs to validation set
  void shift_to_train() const; // shift all of the ptrs to training set
  
  typedef std::vector<int> index_vector;
  index_vector random_order() const;//randomize order of predictor varaiables
  
  double GetBagFraction() const;
  
  long GetValidSize() const;
  long GetTotalInBag() const;

  const bag& GetBag() const
  {
    return dataImpl.afInBag;
  };
  bool GetBagElem(long index) const
  {
	return dataImpl.afInBag[index];
  }
  void FillRemainderOfBag(long offset);
  void SetBagElem(long index, bool value);

private:
  //-------------------
  // Private Variables
  //-------------------
  CDImpl dataImpl;
};

#endif // __dataset_h__


