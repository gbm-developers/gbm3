//------------------------------------------------------------------------------
//
//  File:       gbm.cpp
//
//  Description: sets up distribution for GBM and tree conversion to R.
//
//  History: Greg Ridgeway 2003.
//------------------------------------------------------------------------------

//------------------------------
// Includes
//------------------------------
#include "gbm_setup.h"
#include "distributionFactory.h"
#include <algorithm>
#include <vector>

CDistribution* gbm_setup
(
 const CDataset& data,
 SEXP radMisc,
 const std::string& family,
 int cTrain,
 int& cGroups
 )
{
	cGroups = -1;
	DistributionFactory* factory = DistributionFactory::Get();

	// Checks for pairwise distribution
	// this should be removed later.
	if(0 == family.compare(0, 8, "pairwise"))
	{
		std::size_t offsetToMeasure = family.find("_");
		if(offsetToMeasure == std::string::npos)
		{
			throw GBM::failure("Unable to locate IR metric required for pairwise");
		}

		const char* szIRMeasure = family.c_str() + offsetToMeasure + 1;
		CDistribution* pTemp(factory -> CreateDist("pairwise", radMisc, data, szIRMeasure, cGroups, cTrain));
		return pTemp;
	}else
	{
		CDistribution* pTemp(factory -> CreateDist(family, radMisc, data, "", cGroups, cTrain));
		return pTemp;
	}

}




