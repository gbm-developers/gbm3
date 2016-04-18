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
		adOffset(radOffset), adWeight(radWeight), adX(radX),
		acVarClasses(racVarClasses), alMonotoneVar(ralMonotoneVar),
		aiXOrder(raiXOrder), numOfTrainData(cTrain), numOfFeatures(cFeatures),
		fHasOffset(GBM_FUNC::has_value(adOffset))
	{

		// If you've no offset set to 0
		if(!fHasOffset)
		{
			std::fill(adOffset.begin(), adOffset.begin() + adX.nrow(), 0.0);
		}

		//
		Rcpp::NumericMatrix tempY(radY);
		adY = tempY(Rcpp::_, 1);

		// Set other stuff
		bagFraction = fractionInBag;
		totalInBag = (long) (fractionInBag * cTrain);
		cValid = adX.nrow() - cTrain;
		pointAtTrainSet = true;
		adYPtr = adY.begin();
		adWeightPtr = adWeight.begin();
		adOffsetPtr = adOffset.begin();
		afInBag.assign(cTrain, false);

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

		//-------------------
		// Public Variables
		//-------------------
		// Numeric vectors storing data
		Rcpp::NumericVector adY, adOffset, adWeight;
		Rcpp::NumericMatrix adX;
		Rcpp::IntegerVector acVarClasses, alMonotoneVar, aiXOrder;

		// Ptrs to numeric vectors - these must be mutable
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
//  radOffset      - R object (SEXP) containing the offset applied to each prediction
//  radX           - R object (SEXP) containing the predictor values.
//  raiXOrder      - R object (SEXP) containing the order of predictor values to be
//  used in GBM formula
//  radWeight      - R object (SEXP) containing weights to be used in fitting process
//  racVarClasses  - R object (SEXP) containing the variable classes
//  ralMonotoneVar - R object (SEXP) containing +-1/9 indicating whether the variables are:
//  monotone increasing (+1), decreasing (-1) or arbitrary (0) with the response variable.
//	cTrain		   - int specifiy the number of data points in training set
//
//-----------------------------------
CDataset::CDataset(DataDistParams dataParams)
{
	// Set up the pimpl
	dataImpl = new CDImpl(dataParams.respY, dataParams.offset, dataParams.xValues,
							dataParams.xOrder,	dataParams.varWeight, dataParams.varClasses,
							dataParams.monotoneVar, dataParams.cTrain,
								dataParams.cFeatures, dataParams.dBagFraction);

	// Check for errors on initialization
	if (dataImpl-> adX.ncol() != dataImpl-> alMonotoneVar.size())
	{
		throw GBM::invalid_argument("shape mismatch (monotone does not match data)");
	}

	if (dataImpl-> adX.ncol() != dataImpl-> acVarClasses.size())
	{
		throw GBM::invalid_argument("shape mismatch (var classes does not match data)");
	}

	if (dataImpl->adX.nrow() < int(dataParams.cTrain))
	{
		throw GBM::invalid_argument("your training instances don't make sense");
	}
 };

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
CDataset::~CDataset()
{
	delete dataImpl;
}

//-----------------------------------
// Function: nrow
//
// Returns: number of rows in data
//
// Parameters: none
//
//-----------------------------------
int CDataset::nrow() const
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
int CDataset::ncol() const
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
double* CDataset::y_ptr()
{
	return dataImpl->adYPtr;
}
const double* CDataset::y_ptr() const
{
	return dataImpl->adYPtr;
}

//-----------------------------------
// Function: offset_ptr
//
// Returns: get ptr to offset container
//
// Parameters: none
//
//-----------------------------------
const double* CDataset::offset_ptr(bool require) const
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

double* CDataset::offset_ptr(bool require)
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
double* CDataset::weight_ptr()
{
	return dataImpl->adWeightPtr;
}

const double* CDataset::weight_ptr() const
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
int CDataset::varclass(int ind) const
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
int CDataset::monotone(int ind) const
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
int* CDataset::order_ptr()
{
	return dataImpl->aiXOrder.begin();
}

const int* CDataset::order_ptr() const
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
bool CDataset::has_offset() const
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
double CDataset::x_value(const int row, const int col) const
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
long CDataset::get_trainSize() const
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
long CDataset::get_numFeatures() const
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
void CDataset::shift_to_validation() const
{
	if(dataImpl->pointAtTrainSet)
	{
		dataImpl->adYPtr= dataImpl->shift_ptr_to_validation(dataImpl->adYPtr);
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
void CDataset::shift_to_train() const
{
	if(!(dataImpl->pointAtTrainSet))
	{
		dataImpl->adYPtr = dataImpl->shift_ptr_to_train(dataImpl->adYPtr);
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
index_vector CDataset::random_order() const
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
double CDataset::GetBagFraction() const
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
long CDataset::GetValidSize() const
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
long CDataset::GetTotalInBag() const
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
void CDataset::SetBagElem(long index, bool value)
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
bag CDataset::GetBag()
{
	return dataImpl->afInBag;
}

bag CDataset::GetBag() const
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
bool CDataset::GetBagElem(long index) const
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
void CDataset::FillRemainderOfBag(long offset)
{
	std::fill((dataImpl->afInBag).begin() + offset, (dataImpl->afInBag).end(), false);
}






