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
	      unsigned long whichvar, unsigned long numvar_classes,
	      long monotone);
  
  VarSplitter(const VarSplitter& rhs) :
    initial_sumresiduals_(rhs.initial_sumresiduals_),
    initial_totalweight_(rhs.initial_totalweight_),
    initial_numobs_(rhs.initial_numobs_),
    bestsplit_(rhs.bestsplit_),
    proposedsplit_(rhs.proposedsplit_) {
    splitter_.reset(rhs.splitter_->clone());
  }
  
  VarSplitter& operator=(const VarSplitter& rhs) {
    if (this == &rhs) {
      return *this;
    }
    initial_sumresiduals_ = rhs.initial_sumresiduals_;
    initial_totalweight_ = rhs.initial_totalweight_;
    initial_numobs_ = rhs.initial_numobs_;
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
  double initial_sumresiduals_;
  double initial_totalweight_;
  unsigned long initial_numobs_;

  NodeParams bestsplit_, proposedsplit_;
  std::auto_ptr<generic_splitter_strategy> splitter_;
};
#endif  // VARSPLITTER_H
