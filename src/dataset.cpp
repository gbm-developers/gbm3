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
  dataImpl(dataParams.response, dataParams.offset,
	   dataParams.xvalues,
	   dataParams.xorder, dataParams.variable_weight,
	   dataParams.variable_num_classes,
	   dataParams.variable_monotonicity, dataParams.num_trainrows,
	   dataParams.num_features, dataParams.bagfraction,
	   dataParams.num_trainobservations, dataParams.patientids) {
  
  // Check for errors on initialization
  if (dataImpl.xmatrix_.ncol() != dataImpl.variable_monotonicity_.size())
    {
      throw GBM::InvalidArgument("shape mismatch (monotone does not match data)");
    }

  if (dataImpl.xmatrix_.ncol() != dataImpl.num_variable_classes_.size())
    {
      throw GBM::InvalidArgument("shape mismatch (var classes does not match data)");
    }
  
  if (dataImpl.xmatrix_.nrow() < int(dataParams.num_trainrows))
    {
      throw GBM::InvalidArgument("your training instances don't make sense");
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
CDataset::~CDataset()
{
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
index_vector CDataset::RandomOrder() const
{
  index_vector result(ncol());

  // fill the vector
  for (index_vector::size_type ind=0; ind!=result.size(); ++ind)
    {
      result[ind] = ind;
    }
  
  // and now shuffle
  std::random_shuffle(result.begin(), result.end(), GBM_FUNC::PtrShuffler);
  // and return
  return result;
}







