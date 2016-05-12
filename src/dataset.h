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
#include <memory>
#include <vector>
#include <Rcpp.h>

typedef std::vector<int> bag;

//------------------------------
// Class definition
//------------------------------
class CDataset
{
public:
	//----------------------
	// Public Constructors
	//----------------------
	CDataset(DataDistParams dataParams);

	//---------------------
	// Public destructor
	//---------------------
	virtual ~CDataset();

	//---------------------
	// Public Functions
	//---------------------
	int nrow() const;
	int ncol() const;

	double* y_ptr(long colIndex=0); //get iterator to class labels
	const double* y_ptr(long colIndex=0) const; //const overloaded version

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

	long get_trainSize() const; // get size of training set
	long get_numFeatures() const; // get the number of features in data

	void shift_to_validation() const; // shift all of the ptrs to validation set
	void shift_to_train() const; // shift all of the ptrs to training set

	typedef std::vector<int> index_vector;
	index_vector random_order() const;//randomize order of predictor varaiables
  
	double GetBagFraction() const;

	long GetValidSize() const;
	long GetTotalInBag() const;
	bag GetBag();
	bag GetBag() const;
	bool GetBagElem(long index) const;
	void FillRemainderOfBag(long offset);
	void SetBagElem(long index, bool value);

private:
	//-------------------
	// Private Variables
	//-------------------
	class CDImpl;
	CDImpl* dataImpl;
};


//-----------------------------------
// Class Definition - Private Variable
//-----------------------------------
class CDataset::CDImpl
{
public:
	//----------------------
	// Public Constructors
	//----------------------
	CDImpl(SEXP radY, SEXP radOffset, SEXP radX, SEXP raiXOrder,
		SEXP radWeight, SEXP racVarClasses, SEXP ralMonotoneVar,
		const int cTrain, const int cFeatures, const double fractionInBag):
		adY(radY), adOffset(radOffset), adWeight(radWeight), adX(radX),
		acVarClasses(racVarClasses), alMonotoneVar(ralMonotoneVar),
		aiXOrder(raiXOrder), numOfTrainData(cTrain), numOfFeatures(cFeatures),
		fHasOffset(GBM_FUNC::has_value(adOffset))
	{

		// If you've no offset set to 0
		if(!fHasOffset)
		{
			std::fill(adOffset.begin(), adOffset.begin() + adX.nrow(), 0.0);
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
		inline T* shift_ptr_to_validation(T* x) const
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
		inline T* shift_ptr_to_train(T* x) const
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
		inline void SetUpYPtrs()
		{
			for(long i = 0; i < adY.ncol(); i++)
			{
				yptrs.push_back(adY(Rcpp::_, i).begin());
			}
		}

		//-------------------
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
		bool fHasOffset;
		bool pointAtTrainSet;

		// Bagged  data
		bag afInBag;
		double bagFraction;
		long totalInBag;


};

//-----------------------------------
// Function: nrow
//
// Returns: number of rows in data
//
// Parameters: none
//
//-----------------------------------
inline int CDataset::nrow() const
{
	return dataImpl->adX.nrow();
}

//-----------------------------------
// Function: ncol
//
// Returns: number of columns in data
//
// Parameters: none
//
//-----------------------------------
inline int CDataset::ncol() const
{
	return dataImpl->adX.ncol();
}

//-----------------------------------
// Function: y_ptr
//
// Returns: iterator to response variable container
//
// Parameters: none
//
//-----------------------------------
inline double* CDataset::y_ptr(long colIndex)
{
	return dataImpl->yptrs[colIndex];
}
inline const double* CDataset::y_ptr(long colIndex) const
{
	return dataImpl->yptrs[colIndex];
}

//-----------------------------------
// Function: offset_ptr
//
// Returns: get ptr to offset container
//
// Parameters: none
//
//-----------------------------------
inline const double* CDataset::offset_ptr(bool require) const
{
	if (has_offset())
	{
	  return dataImpl->adOffsetPtr;
	}
	else
	{
	  if (require)
	  {
		throw GBM::failure("You require a genuine offset, and don't have one.");
	  }
	  else
	  {
		return 0;
	  }
	}
}

inline double* CDataset::offset_ptr(bool require)
{
	return const_cast<double*>(static_cast<const CDataset*>(this)->offset_ptr(require));
}

//-----------------------------------
// Function: weight_ptr
//
// Returns: get ptr to weights conta	unsigned long GetValidSize();iner
//
// Parameters: none
//
//-----------------------------------
inline double* CDataset::weight_ptr()
{
	return dataImpl->adWeightPtr;
}

inline const double* CDataset::weight_ptr() const
{
	return dataImpl->adWeightPtr;
}

//-----------------------------------
// Function: varclass
//
// Returns: variable class - an int
//
// Parameters: ind - integer specifying index of data-point whose class to look-up.
//
//-----------------------------------
inline int CDataset::varclass(int ind) const
{
	return dataImpl->acVarClasses[ind];
}

//-----------------------------------
// Function: monotone
//
// Returns: int indicating relationship between predictor and response.
//
// Parameters: ind - integer specifying index of data-point
//
//-----------------------------------
inline int CDataset::monotone(int ind) const
{
	return dataImpl->alMonotoneVar[ind];
}

//-----------------------------------
// Function: order_ptr
//
// Returns: iterator to order of predictors container.
//
// Parameters: none
//
//-----------------------------------
inline int* CDataset::order_ptr()
{
	return dataImpl->aiXOrder.begin();
}

inline const int* CDataset::order_ptr() const
{
	return dataImpl->aiXOrder.begin();
}

//-----------------------------------
// Function: has_offset
//
// Returns: bool indicating if the function has an offset.
//
// Parameters: none
//
//-----------------------------------
inline bool CDataset::has_offset() const
{
	return dataImpl->fHasOffset;
}

//-----------------------------------
// Function: x_value
//
// Returns: double - the predictor value desired.
//
// Parameters:
//  row - int representing the row of the data-point in the array
//  col - int representing the row of the data-point in the array
//
//-----------------------------------
inline double CDataset::x_value(const int row, const int col) const
{
	return dataImpl->adX(row, col);
}

//-----------------------------------
// Function: get_trainSize
//
// Returns: long - the number of training instances
//
// Parameters: none
//
//-----------------------------------
inline long CDataset::get_trainSize() const
{
	return dataImpl->numOfTrainData;
}

//-----------------------------------
// Function: get_numFeatures
//
// Returns: long - the number of features
//
// Parameters: none
//
//-----------------------------------
inline long CDataset::get_numFeatures() const
{
	return dataImpl->numOfFeatures;
}

//-----------------------------------
// Function: shift_to_validation
//
// Returns:  shifts the data to the validation set.
//
// Parameters: none
//
//-----------------------------------
inline void CDataset::shift_to_validation() const
{
	if(dataImpl->pointAtTrainSet)
	{
		for(int i = 0; i < dataImpl->yptrs.size(); i++)
		{
			dataImpl->yptrs[i] = dataImpl->shift_ptr_to_validation(dataImpl->yptrs[i]);
		}
		dataImpl->adOffsetPtr = dataImpl->shift_ptr_to_validation(dataImpl->adOffsetPtr);
		dataImpl->adWeightPtr = dataImpl->shift_ptr_to_validation(dataImpl->adWeightPtr);
		dataImpl->pointAtTrainSet = false;
	}
	else
	{
		throw GBM::invalid_argument("Data is already the validation set.");
	}
}

//-----------------------------------
// Function: shift_to_train
//
// Returns:  shifts the data to the training set.
//
// Parameters: none
//
//-----------------------------------
inline void CDataset::shift_to_train() const
{
	if(!(dataImpl->pointAtTrainSet))
	{
		for(int i = 0; i < dataImpl->yptrs.size(); i++)
		{
			dataImpl->yptrs[i] = dataImpl->shift_ptr_to_train(dataImpl->yptrs[i]);
		}
		dataImpl->adOffsetPtr = dataImpl->shift_ptr_to_train(dataImpl->adOffsetPtr);
		dataImpl->adWeightPtr = dataImpl->shift_ptr_to_train(dataImpl->adWeightPtr);
		dataImpl->pointAtTrainSet = true;
	}
	else
	{
		throw GBM::invalid_argument("Data is already the training set.");
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
typedef std::vector<int> index_vector;
inline index_vector CDataset::random_order() const
{
	index_vector result(ncol());

	// fill the vector
	for (index_vector::size_type ind=0; ind!=result.size(); ++ind)
	{
		result[ind] = ind;
	}

	// and now shuffle
	std::random_shuffle(result.begin(), result.end(), GBM_FUNC::ptrShuffler);
	// and return
	return result;
}

//-----------------------------------
// Function: GetBagFraction
//
// Returns: double
//
// Description: get fraction of data in bag
//
// Parameters: none
//
//-----------------------------------
inline double CDataset::GetBagFraction() const
{
	return dataImpl->bagFraction;
}

//-----------------------------------
// Function: GetValidSize
//
// Returns: long
//
// Description: get size of validation set.
//
// Parameters: none
//
//-----------------------------------
inline long CDataset::GetValidSize() const
{
	return dataImpl->cValid;
}

//-----------------------------------
// Function: GetTotalInBag
//
// Returns: long
//
// Description: get total amount of data in bag
//
// Parameters: none
//
//-----------------------------------
inline long CDataset::GetTotalInBag() const
{
	return dataImpl->totalInBag;
}

//-----------------------------------
// Function: SetBagElem
//
// Returns: none
//
// Description: setter for bag elements
//
// Parameters: long - index of element to set
//    bool - value to set index bag to.
//
//-----------------------------------
inline void CDataset::SetBagElem(long index, bool value)
{
	dataImpl->afInBag[index] = value;
}

//-----------------------------------
// Function: GetBag
//
// Returns: bag
//
// Description: getter for bag
//
// Parameters: none
//
//-----------------------------------
inline bag CDataset::GetBag()
{
	return dataImpl->afInBag;
}

inline bag CDataset::GetBag() const
{
	return dataImpl->afInBag;
}

//-----------------------------------
// Function: GetBagElem
//
// Returns: const bool
//
// Description: getter for bag element
//
// Parameters: long - index of element to get
//
//-----------------------------------
inline bool CDataset::GetBagElem(long index) const
{
	return dataImpl->afInBag[index];
}

//-----------------------------------
// Function: FillRemainderOfBag
//
// Returns: void
//
// Description: once data put in bag set all other elems to False
//
// Parameters: long - sets where to start filling the bag
//
//-----------------------------------
inline void CDataset::FillRemainderOfBag(long offset)
{
	std::fill((dataImpl->afInBag).begin() + offset, (dataImpl->afInBag).end(), false);
}

#endif // __dataset_h__


