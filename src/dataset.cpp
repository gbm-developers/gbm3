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
CDataset::CDataset(const DataDistParams& dataParams) :
  dataImpl(dataParams.respY, dataParams.offset,
	   dataParams.xValues,
	   dataParams.xOrder, dataParams.varWeight,
	   dataParams.varClasses,
	   dataParams.monotoneVar, dataParams.cTrain,
	   dataParams.cFeatures, dataParams.dBagFraction) {
  
  // Check for errors on initialization
  if (dataImpl.adX.ncol() != dataImpl.alMonotoneVar.size())
    {
      throw GBM::invalid_argument("shape mismatch (monotone does not match data)");
    }

  if (dataImpl.adX.ncol() != dataImpl.acVarClasses.size())
    {
      throw GBM::invalid_argument("shape mismatch (var classes does not match data)");
    }
  
  if (dataImpl.adX.nrow() < int(dataParams.cTrain))
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
	return dataImpl.numOfTrainData;
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
	return dataImpl.numOfFeatures;
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
  dataImpl.shift_to_validation();
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
  dataImpl.shift_to_train();
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
	return dataImpl.bagFraction;
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
  return dataImpl.cValid;
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
  return dataImpl.totalInBag;
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
	dataImpl.afInBag[index] = value;
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
	return dataImpl.afInBag;
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
	std::fill((dataImpl.afInBag).begin() + offset, (dataImpl.afInBag).end(), false);
}






