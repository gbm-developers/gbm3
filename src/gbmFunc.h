//////////////////////////////////////////////
//
// File: gbmFunc.h
//
// Author: James Hickey
//
// Description: Functions that are accessible to all of GBM.
//
//////////////////////////////////////////////
#include <Rcpp.h>

namespace GBM_FUNC
{
	int numGroups(const double* adMisc, int cTrain);
	bool has_value(const Rcpp::NumericVector& vec);
	std::ptrdiff_t ptrShuffler(std::ptrdiff_t n);
}
