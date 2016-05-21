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

void CDistribution::bagIt(CDataset& data, std::multimap<int, int>& patIdToRow) {

	unsigned long i = 0;
	unsigned long cBagged = 0;

	pair<multimap<int, int>::iterator, multimap<int, int>::iterator> keyRange;
	multimap<int, int>::iterator patIt, rowIt;

	// Bag via patient id  - loop over patients
	for(patIt = patIdToRow.begin();
			patIt != patIdToRow.end();
			patIt = rowIt)
	{

		// Check if we've filled the bag or have left the training set
		if((i >= data.GetNumPatientsInTraining()) || (cBagged >= data.GetTotalInBag())) break;

		keyRange = patIdToRow.equal_range(patIt->first);

		// Check if that patient should be bagged - bag corresponding rows
		if(unif_rand() * (data.GetNumPatientsInTraining()-i) < data.GetTotalInBag() - cBagged)
		{
			cBagged++;
			for(rowIt = keyRange.first; rowIt != keyRange.second; ++rowIt)
			{
				data.SetBagElem((*rowIt).second);
			}
		}
		else
		{
			rowIt = keyRange.second;
		}

		// Increment patient number
		i += 1;
	}
}


