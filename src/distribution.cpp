//  GBM by Greg Ridgeway  Copyright (C) 2003

#include "distribution.h"


CDistribution::CDistribution(SEXP radMisc, const CDataset& data)
: adMisc(radMisc), distHasMisc(has_value(adMisc))
{
	pData = &data;
}

CDistribution::~CDistribution()
{

}

bool CDistribution::has_misc() const
{
   return distHasMisc;
}

const double* CDistribution::misc_ptr(bool require) const
{
	 if (has_misc())
	 {
	   return adMisc.begin();
	 }
	 else
	 {
	   if (require)
	   {
		 throw GBM::failure("You require genuine misc, and don't have it.");
	   }
	   else
	   {
		 return 0;
	   }
	 }
}

double* CDistribution::misc_ptr(bool require)
{
 return const_cast<double*>(static_cast<const CDistribution*>(this)->misc_ptr(require));
}

const CDataset* CDistribution::data_ptr() const
{
	return pData;
}




