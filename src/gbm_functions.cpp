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
int gbm_functions::NumGroups(const double* kMisc, int num_training_rows)
{
	if (num_training_rows <= 0)
	{
		return 0;
	}
	double lastgroup = kMisc[0];
	int groups = 1;

	for(int i=1; i<num_training_rows; i++)
	{
		const double dGroup = kMisc[i];
		if (dGroup != lastgroup)
		{
			lastgroup = dGroup;
			groups++;
		}
	}
	return groups;
}

// Function that checks that a R vector is not NA.
bool gbm_functions::has_value(const Rcpp::NumericVector& kVec)
{
    return !( (kVec.size() == 1) && (ISNA(kVec[0])));
}

// Function that shuffles an array.
std::ptrdiff_t gbm_functions::PtrShuffler(std::ptrdiff_t n)
{
    return n * unif_rand();
}


