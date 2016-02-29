//////////////////////////////////////////////
//
// File: numGroups.cpp
//
// Author: James Hickey
//
// Description: Function that counts the number of distinct
//				groups in the input data. Used only with pairwise distributions.
//
//////////////////////////////////////////////
#include "numGroups.h"

int GBM_Func::numGroups(const double* adMisc, int cTrain)
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

