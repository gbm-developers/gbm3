//------------------------------------------------------------------------------
//  File:    vec_varsplitters.h
//
//	Description: header file for class specifying the vector of var
// splitters
//
//  Created on: 10 Jun 2016
//
//  Author: jhickey
//------------------------------------------------------------------------------

#ifndef VEC_VARSPLITTERS_H_
#define VEC_VARSPLITTERS_H_

//------------------------------
// Includes
//------------------------------
#include "varsplitter.h"
#include "vec_nodeparams.h"
#include <Rcpp.h>

//------------------------------
// Class Definition
//------------------------------
class VecVarSplitters {
 public:
  //----------------------
  // Public Constructors
  //----------------------
  VecVarSplitters(unsigned long memory_space,
                  const std::vector<CNode*>& kVecTermNodePtrs,
                  const unsigned long kMinNumNodeObs, unsigned long bias,
                  const unsigned long kVarNum, const unsigned long kVarClasses,
                  long monotone) {
    varsplitters_.reserve(memory_space);
    for (unsigned long node_num = 0; node_num < memory_space; node_num++) {
      varsplitters_.push_back(VarSplitter(*kVecTermNodePtrs[node_num],
                                          kMinNumNodeObs, bias, kVarNum,
                                          kVarClasses, monotone));
    }
  }

  VecVarSplitters(const VecVarSplitters& rhs)
      : varsplitters_(rhs.varsplitters_) {}

  //---------------------
  // Public destructor
  //---------------------
  ~VecVarSplitters(){};

  //---------------------
  // Public Functions
  //---------------------
  unsigned long size() const { return varsplitters_.size(); };
  VarSplitter& operator[](unsigned long node_num) {
    return varsplitters_[node_num];
  };

  VecNodeParams proposal() {
    VecNodeParams proposed_splits;
    proposed_splits.reserve(varsplitters_.size());

    for (unsigned long node_num = 0; node_num < varsplitters_.size();
         node_num++) {
      proposed_splits.push_back(varsplitters_[node_num].best_split());
    }
    return proposed_splits;
  }

 private:
  //---------------------
  // Private Variables
  //---------------------
  std::vector<VarSplitter> varsplitters_;
};
#endif /* VEC_VARSPLITTERS_H_ */
