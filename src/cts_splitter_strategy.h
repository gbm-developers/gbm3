#ifndef CTS_SPLITTER_STRATEGY_H
#define CTS_SPLITTER_STRATEGY_H

#include "generic_splitter_strategy.h"

class cts_splitter_strategy : public generic_splitter_strategy {
 public:
 cts_splitter_strategy(unsigned long min_num_node_obs,
		       long monotonicity)
   : last_xvalue_(-HUGE_VAL), min_num_node_obs_(min_num_node_obs),
    monotonicity_(monotonicity) {};

  generic_splitter_strategy* clone() const {
    return new cts_splitter_strategy(*this);
  }

  void incorporate_obs(NodeParams& bestsplit, NodeParams& proposedsplit,
                       double xval, double residval, double weight) {
    if (ISNA(xval)) {
      proposedsplit.UpdateMissingNode(weight * residval, weight);
      return;
    }

    if (last_xvalue_ > xval) {
      throw gbm_exception::Failure(
          "Observations are not in order. gbm() was unable to build an index "
          "for the design matrix. Could be a bug in gbm or an unusual data "
          "type in data.");
    }

    // Evaluate the current split
    // the newest observation is still in the right child
    proposedsplit.set_split_value(0.5 * (last_xvalue_ + xval));

    if ((last_xvalue_ != xval) &&
        proposedsplit.has_min_num_obs(min_num_node_obs_) &&
        proposedsplit.split_is_correct_monotonicity(monotonicity_)) {
      proposedsplit.NodeGradResiduals();
      if (proposedsplit.get_improvement() > bestsplit.get_improvement()) {
        bestsplit = proposedsplit;
      }
    }

    // now move the new observation to the left
    // if another observation arrives we will evaluate this
    proposedsplit.UpdateLeftNode(weight * residval, weight);
    last_xvalue_ = xval;
  };

  void wrap_up(NodeParams& bestsplit, NodeParams& proposedsplit) {
    return;
  };

 private:
  double last_xvalue_;
  unsigned long min_num_node_obs_;
  long monotonicity_;
};
#endif
