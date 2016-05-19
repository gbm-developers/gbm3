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

void CDistribution::bagIt(CDataset& data) {

   for(unsigned long i=0, cBagged = 0;
      (i<data.get_trainSize()) && (cBagged < data.GetTotalInBag());
      i++) {
    if ((unif_rand() * (data.get_trainSize() - i)) <
	(data.GetTotalInBag() - cBagged)) {
      data.SetBagElem(i);
      cBagged++;
    }
  }

}

