#ifndef CATEGORICAL_SPLITTER_STRATEGY_H
#define CATEGORICAL_SPLITTER_STRATEGY_H

#include "generic_splitter_strategy.h"

class categorical_splitter_strategy : public generic_splitter_strategy {
 public:
  categorical_splitter_strategy(unsigned long min_num_node_obs,
                                unsigned long size)
      : min_num_node_obs_(min_num_node_obs), group_(size) {};

  generic_splitter_strategy* clone() const {
    return new categorical_splitter_strategy(*this);
  }

  void incorporate_obs(NodeParams& bestsplit, NodeParams& proposedsplit,
                       double xval, double residval, double weight) {
    if (ISNA(xval)) {
      proposedsplit.UpdateMissingNode(weight * residval, weight);
      return;
    }

    unsigned long cat = xval;
    group_[cat].increment(weight * residval, weight, 1);
  }

  void wrap_up(NodeParams& bestsplit, NodeParams& proposedsplit) {
    std::vector<std::pair<double, int> > groupMeanAndCat(group_.size());
    unsigned long num_finite_means = 0;

    // sort the groups

    for (std::size_t ind = 0; ind != groupMeanAndCat.size(); ++ind) {
      groupMeanAndCat[ind].second = ind;
      if (group_[ind].get_totalweight() > 0) {
        groupMeanAndCat[ind].first = group_[ind].prediction();
	num_finite_means++;
      } else {
        groupMeanAndCat[ind].first = HUGE_VAL;
      }
    }

    std::sort(groupMeanAndCat.begin(), groupMeanAndCat.end());

    for (std::size_t ind = 0;
         (num_finite_means > 1) && (1 + ind < num_finite_means); ind++) {
      proposedsplit.set_split_value(ind);

      const NodeDef& weights = group_[groupMeanAndCat[ind].second];
      proposedsplit.UpdateLeftNode(weights.get_weightresid(),
                                   weights.get_totalweight(),
                                   weights.get_num_obs());
      proposedsplit.NodeGradResiduals();

      if (proposedsplit.has_min_num_obs(min_num_node_obs_) &&
          (proposedsplit.get_improvement() > bestsplit.get_improvement())) {
        bestsplit = proposedsplit;
      }
    }

    // let's not copy this lots
    if (num_finite_means > 1) {
      bestsplit.SetBestCategory(groupMeanAndCat);
    }
  };

 private:
  unsigned long min_num_node_obs_;
  std::vector<NodeDef> group_;
};

#endif
