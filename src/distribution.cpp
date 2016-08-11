//  GBM by Greg Ridgeway  Copyright (C) 2003

#include "distribution.h"

CDistribution::CDistribution() : parallel_(), num_groups_(-1) {}

CDistribution::CDistribution(const parallel_details& parallel)
    : parallel_(parallel), num_groups_(-1) {}

CDistribution::~CDistribution() {}

void CDistribution::BagData(const CDataset& kData, Bag& bag) {
  unsigned long i = 0;
  unsigned long numbagged = 0;

  pair<multimap<int, int>::iterator, multimap<int, int>::iterator> keyrange;
  multimap<int, int>::iterator obs_it, row_it;

  // Bag via patient id  - loop over observations
  for (obs_it = obsid_to_row_.begin(); obs_it != obsid_to_row_.end();
       obs_it = obsid_to_row_.upper_bound(obs_it->first)) {

	// Check if we've filled the bag or have left the training set
    // Works as long as ids are sequential
    if ((i >= kData.get_num_observations_in_training()) ||
        (numbagged >= bag.get_total_in_bag()))
      break;

    keyrange = obsid_to_row_.equal_range(obs_it->first);

    // Check if that observation should be bagged - bag corresponding rows
    if (unif_rand() * (kData.get_num_observations_in_training() - i) <
        bag.get_total_in_bag() - numbagged) {
      numbagged++;
      for (row_it = keyrange.first; row_it != keyrange.second; ++row_it) {
    	  bag.set_element((*row_it).second);
      }
    }

    // Increment observation number
    i += 1;
  }
}
