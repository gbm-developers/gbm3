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
  VecVarSplitters(unsigned long memory_space) { varsplitters_.reserve(memory_space); };
  VecVarSplitters(unsigned long memory_space, const VarSplitter& varsplitter) :
	  varsplitters_(memory_space, varsplitter) {};

  //---------------------
  // Public destructor
  //---------------------
  ~VecVarSplitters() {};

  //---------------------
  // Public Functions
  //---------------------
  void push_back(const VarSplitter& varsplitter) {
	  varsplitters_.push_back(varsplitter);
  };
  unsigned long size() const { return varsplitters_.size(); };
  VarSplitter& operator[](unsigned long node_num) { return varsplitters_[node_num]; };
  VecVarSplitters& operator+=(const VecVarSplitters& rhs) {
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
