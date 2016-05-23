//////////////////////////////////////////////
//
// File: gbmFunc.cpp
//
// Author: James Hickey
//
// Description: Functions that are accessible to all of GBM.
//
//////////////////////////////////////////////
#include "gbm_functions.h"

//  Function that counts the number of distinct groups in the input data.
//  Used only with the pairwise distribution.
int GBM_FUNC::numGroups(const double* adMisc, int cTrain)
{
	if (cTrain <= 0)
	{
		return 0;
	}
	double dLastGroup = adMisc[0];
	int cGroups = 1;

	for(int i=1; i<cTrain; i++)
	{
		const double dGroup = adMisc[i];
		if (dGroup != dLastGroup)
		{
			dLastGroup = dGroup;
			cGroups++;
		}
	}
	return cGroups;
}

// Function that checks that a R vector is not NA.
bool GBM_FUNC::has_value(const Rcpp::NumericVector& x)
{
    return !( (x.size() == 1) && (ISNA(x[0])));
}

// Function that shuffles an array.
std::ptrdiff_t GBM_FUNC::ptrShuffler(std::ptrdiff_t n)
{
    return n * unif_rand();
}


