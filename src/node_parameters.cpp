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
#include "gbmexcept.h"

//---------------------
// Public Functions
//---------------------
NodeParams::~NodeParams() {}

void NodeParams::ResetSplitProperties(double weightedresiduals,
                                      double trainingweight,
                                      unsigned long numobs, double splitvalue,
                                      unsigned long variableclasses,
                                      unsigned long splitvar) {
  right_ = NodeDef(weightedresiduals, trainingweight, numobs);

  left_.clear();
  missing_.clear();

  split_var_ = splitvar;
  split_value_ = splitvalue;
  improvement_ = 0.0;
  split_class_ = variableclasses;
}
