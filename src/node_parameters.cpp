//------------------------------------------------------------------------------
//
//  File:       nodeParameters.h
//
//  Description: class implementing the parameters in each node
//
//------------------------------------------------------------------------------
//---------------------
// Includes
//---------------------
#include "node_parameters.h"

//---------------------
// Public Functions
//---------------------
NodeParams::~NodeParams() {}

void NodeParams::ResetSplitProperties(double weightedresiduals,
                                      double trainingweight,
                                      unsigned long numobs,
                                      unsigned long variableclasses,
                                      unsigned long splitvar) {
  right_ = NodeDef(weightedresiduals, trainingweight, numobs);

  left_.clear();
  missing_.clear();

  split_var_ = splitvar;
  split_value_ = -HUGE_VAL;
  improvement_ = 0.0;
  split_class_ = variableclasses;
}
