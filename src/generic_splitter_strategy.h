#ifndef GENERIC_SPLITTER_STRATEGY_H
#define GENERIC_SPLITTER_STRATEGY_H

#include "node_parameters.h"

// and here's the class

class generic_splitter_strategy {
 public:
  //------------------
  // Public destructor
  //------------------
  virtual ~generic_splitter_strategy(){};
  virtual generic_splitter_strategy* clone() const = 0;
  virtual void incorporate_obs(NodeParams& bestsplit, NodeParams& proposedsplit,
                               double xval, double residval, double weight) = 0;
  virtual void wrap_up(NodeParams& bestsplit, NodeParams& proposedsplit) = 0;
};

#endif
