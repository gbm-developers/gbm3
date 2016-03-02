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
#include "gbm.h"
#include "distributionFactory.h"
#include <algorithm>
#include <vector>

std::auto_ptr<CDistribution> gbm_setup
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
		std::auto_ptr<CDistribution> pDist(factory -> CreateDist("pairwise", radMisc, data, szIRMeasure, cGroups, cTrain));
		return pDist;
	}else{
		std::auto_ptr<CDistribution> pDist(factory -> CreateDist(family, radMisc, data, "", cGroups, cTrain));
		return pDist;
	}

}


void gbm_transfer_to_R
(
 CGBM *pGBM,
 VEC_VEC_CATEGORIES &vecSplitCodes,
 int *aiSplitVar,
 double *adSplitPoint,
 int *aiLeftNode,
 int *aiRightNode,
 int *aiMissingNode,
 double *adErrorReduction,
 double *adWeight,
 double *adPred,
 int cCatSplitsOld
 )
{
    pGBM->TransferTreeToRList(aiSplitVar,
			      adSplitPoint,
			      aiLeftNode,
			      aiRightNode,
			      aiMissingNode,
			      adErrorReduction,
			      adWeight,
			      adPred,
			      vecSplitCodes,
			      cCatSplitsOld);
}




