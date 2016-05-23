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
  right = NodeDef(weightedResiduals, trainingWeight, numObs);

  left.clear();
  missing.clear();
  
  SplitVar = splitVar;
  SplitValue = splitValue;
  ImprovedResiduals = 0.0;
  SplitClass = variableClasses;
}
