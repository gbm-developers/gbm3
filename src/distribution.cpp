//  GBM by Greg Ridgeway  Copyright (C) 2003

#include "distribution.h"


CDistribution::CDistribution()
{
	cGroups = -1;
}

CDistribution::~CDistribution()
{

}

int CDistribution::GetNumGroups() const
{
	return cGroups;
}

void CDistribution::SetNumGroups(int GroupVal)
{
	cGroups = GroupVal;
}




