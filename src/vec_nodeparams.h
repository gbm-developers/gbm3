//------------------------------------------------------------------------------
//  File:    vec_nodeparams.h
//
//	Description: header file for class specifying the vector of node params
//
//  Created on: 14 Jun 2016
//
//  Author: jhickey
//------------------------------------------------------------------------------

#ifndef VEC_NODEPARAMS_H_
#define VEC_NODEPARAMS_H_

//------------------------------
// Includes
//------------------------------
#include "node_parameters.h"

//------------------------------
// Class Definition
//------------------------------
class VecNodeParams {
 public:
  //----------------------
  // Public Constructors
  //----------------------
  VecNodeParams(){};
  VecNodeParams(unsigned long memory_space) : nodeparams_(memory_space) {};
  VecNodeParams(const VecNodeParams& rhs) : nodeparams_(rhs.nodeparams_) {};
  //---------------------
  // Public destructor
  //---------------------
  ~VecNodeParams(){};

  //---------------------
  // Public Functions
  //---------------------
  unsigned long size() const { return nodeparams_.size(); };
  void reserve(unsigned long memory_space) {
    nodeparams_.reserve(memory_space);
  };
  void push_back(const NodeParams& nodeparams) {
    nodeparams_.push_back(nodeparams);
  };
  NodeParams& operator[](unsigned long node_num) {
    return nodeparams_[node_num];
  };
  VecNodeParams& operator+=(const VecNodeParams& rhs) {
    if (rhs.size() > size()) {
      throw gbm_exception::Failure(
          "VecNodeParams do not"
          " have compatible sizes");
    }
    for (unsigned long node_num = 0; node_num < rhs.size(); node_num++) {
      nodeparams_[node_num] += rhs.nodeparams_[node_num];
    }
    return *this;
  }

 private:
  //---------------------
  // Private Variables
  //---------------------
  std::vector<NodeParams> nodeparams_;
};
#endif /* VEC_NODEPARAMS_H_ */
