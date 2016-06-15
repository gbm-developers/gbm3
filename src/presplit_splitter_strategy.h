#ifndef PRESPLIT_SPLITTER_STRATEGY_H
#define PRESPLIT_SPLITTER_STRATEGY_H

#include "generic_splitter_strategy.h"

class presplit_splitter_strategy : public generic_splitter_strategy {
public:
  generic_splitter_strategy* clone() const {
    return new presplit_splitter_strategy();
  };
  void incorporate_obs(NodeParams& bestsplit,
		       NodeParams& proposedsplit,
		       double xval, double residval, double weight) { return; };
  void wrap_up(NodeParams& bestsplit,
	       NodeParams& proposedsplit) { return; };
};

#endif
