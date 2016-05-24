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
NodeParams::~NodeParams()
{

}

void NodeParams::ResetSplitProperties(double weightedResiduals,
				      double trainingWeight,
				      unsigned long numObs,
				      double splitValue,
				      unsigned long variableClasses,
				      unsigned long splitVar)
{
  right_ = NodeDef(weightedResiduals, trainingWeight, numObs);

  left_.clear();
  missing_.clear();
  
  split_var_ = splitVar;
  split_value_ = splitValue;
  improvement_ = 0.0;
  split_class_ = variableClasses;
}
