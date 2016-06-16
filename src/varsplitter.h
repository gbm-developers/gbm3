//------------------------------------------------------------------------------
//
//  File:       varsplitter.h
//
//  Description: header for class that splits a node on a particular variable.
//
//------------------------------------------------------------------------------

#ifndef VARSPLITTER_H
#define VARSPLITTER_H

//------------------------------
// Includes
//------------------------------
#include "node.h"
#include "node_parameters.h"
#include "generic_splitter_strategy.h"
#include <memory>

//------------------------------
// Class Definition
//------------------------------
class VarSplitter {
 public:
  //----------------------
  // Public Constructors
  //----------------------
  VarSplitter(CNode& nodetosplit, unsigned long min_num_node_obs,
              unsigned long bias, unsigned long whichvar,
              unsigned long numvar_classes, long monotone);

  VarSplitter(const VarSplitter& rhs)
      : initial_(rhs.initial_),
        bestsplit_(rhs.bestsplit_),
        proposedsplit_(rhs.proposedsplit_) {
    splitter_.reset(rhs.splitter_->clone());
  }

  VarSplitter& operator=(const VarSplitter& rhs) {
    if (this == &rhs) {
      return *this;
    }
    initial_ = rhs.initial_;
    bestsplit_ = rhs.bestsplit_;
    proposedsplit_ = rhs.proposedsplit_;
    splitter_.reset(rhs.splitter_->clone());
    return *this;
  };

  //---------------------
  // Public Functions
  //---------------------

  void IncorporateObs(double xval, double residval, double weight);
  void WrapUpCurrentVariable();

  const NodeParams& best_split() const { return bestsplit_; };

 private:
  //---------------------
  // Private Variables
  //---------------------
  NodeDef initial_;

  NodeParams bestsplit_, proposedsplit_;
  std::auto_ptr<generic_splitter_strategy> splitter_;
};
#endif  // VARSPLITTER_H
