//------------------------------------------------------------------------------
//  File:    vec_varsplitters.h
//
//	Description: header file for class specifying the vector of
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
#include <Rcpp.h>

//------------------------------
// Class Definition
//------------------------------
class VecVarSplitters {
 public:
  //----------------------
  // Public Constructors
  //----------------------
  VecVarSplitters() {};
  VecVarSplitters(unsigned long memory_space) { varsplitters_.reserve(memory_space); };
  VecVarSplitters(unsigned long memory_space, const VarSplitter& varsplitter) :
	  varsplitters_(memory_space, varsplitter) {};

  VecVarSplitters(unsigned long memory_space,
		  	  	  const std::vector<CNode*>& kVecTermNodePtrs,
		  	  	  const unsigned long kMinNumNodeObs,
		  	  	  const unsigned long kVarNum, const unsigned long kVarClasses) {
	  varsplitters_.reserve(memory_space);
	  for (unsigned long node_num = 0; node_num < memory_space;
	          node_num++) {
	     	varsplitters_.push_back(VarSplitter(*kVecTermNodePtrs[node_num],
	     					kMinNumNodeObs, kVarNum, kVarClasses));
	   }
  }

  VecVarSplitters(const VecVarSplitters& rhs) {
	  varsplitters_.reserve(rhs.size());
	  for(unsigned long node_num = 0; node_num < rhs.size(); node_num++) {
		  varsplitters_.push_back(rhs.varsplitters_[node_num]);
	  }
  }

  //---------------------
  // Public destructor
  //---------------------
  ~VecVarSplitters() {};

  //---------------------
  // Public Functions
  //---------------------
  void reserve(const unsigned long memory_space) {
	  varsplitters_.reserve(memory_space);
  }
  void push_back(const VarSplitter& varsplitter) {
	  varsplitters_.push_back(varsplitter);
  };
  unsigned long size() const { return varsplitters_.size(); };
  VarSplitter& operator[](unsigned long node_num) { return varsplitters_[node_num]; };
  VecVarSplitters& operator+=(const VecVarSplitters& rhs) {
	  if(rhs.size() > size()) {
		  throw gbm_exception::Failure("VecVarSplitters do not"
				  " have compatibale sizes");
	  }
	  for(unsigned long node_num = 0; node_num < rhs.size(); node_num++) {
		  varsplitters_[node_num] += rhs.varsplitters_[node_num];
	  }
	  return *this;
  }


 private:
  //---------------------
  // Private Variables
  //---------------------
  std::vector<VarSplitter> varsplitters_;
};
#endif /* VEC_VARSPLITTERS_H_ */
