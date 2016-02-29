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
//
//------------------------------------------------------------------------------

#ifndef __dataset_h__
#define __dataset_h__

//------------------------------
// Includes
//------------------------------
#include "buildinfo.h"
#include "gbmexcept.h"
#include<memory>
#include <vector>
#include <Rcpp.h>

typedef std::vector<int> bag;

namespace {
  inline bool has_value(const Rcpp::NumericVector& x) {
    return !( (x.size() == 1) && (ISNA(x[0])));
  }

  std::ptrdiff_t shuffler(std::ptrdiff_t n) {
    return n * unif_rand();
  }
}
  
//------------------------------
// Class definition
//------------------------------
class CDataset
{
public:
	//----------------------
	// Public Constructors
	//----------------------
	CDataset(SEXP radY, SEXP radOffset, SEXP radX, SEXP raiXOrder,
			SEXP radWeight, SEXP racVarClasses, SEXP ralMonotoneVar, const int cTrain);

	//---------------------
	// Public destructor
	//---------------------
	virtual ~CDataset();

	//---------------------
	// Public Functions
	//---------------------
	int nrow() const;
	int ncol() const;

	double* y_ptr(); //get iterator to class labels
	const double* y_ptr() const; //const overloaded version

	double* offset_ptr(bool require=true); //get iterator to fit offset
	const double* offset_ptr(bool require=true) const; //const overloaded version

	double* weight_ptr(); //get iterator to weights on each data-point
	const double* weight_ptr() const; //const overloaded version

	int varclass(int ind) const;
	int monotone(int ind) const;

	int* order_ptr();//get iterator to ordering of predictor variables
	const int* order_ptr() const;//const overloaded version

	bool has_offset() const;
	double x_value(const int row, const int col) const; // retrieve predictor value
	unsigned long get_trainSize() const; // get size of training set

	void shift_to_validation() const; // shift all of the ptrs to validation set
	void shift_to_train() const; // shift all of the ptrs to training set

	typedef std::vector<int> index_vector;
	index_vector random_order() const;//randomize order of predictor varaiables
  
private:
	//-------------------
	// Private Variables
	//-------------------
	class CDImpl;
	std::auto_ptr<CDImpl> dataImpl;
};

#endif // __dataset_h__


