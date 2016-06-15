//-----------------------------------
//
// File: varsplitter.cpp
//
// Description: class that implements the splitting of a node on a variable.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "varsplitter.h"
#include "splitter_strategies.h"

//---------------------
// Public Functions
//---------------------

VarSplitter::VarSplitter(CNode& nodetosplit,
			 unsigned long min_num_node_obs,
			 unsigned long whichvar,
			 unsigned long numvar_classes,
			 long monotone)
  : initial_(nodetosplit.as_node_def()),
    bestsplit_(initial_),
    proposedsplit_(initial_, numvar_classes, whichvar) {

  if (nodetosplit.is_split_determined()) {
    splitter_.reset(new presplit_splitter_strategy());
  } else if (proposedsplit_.split_class() == 0) {
    // continuous variable
    splitter_.reset(new cts_splitter_strategy(min_num_node_obs,
					      monotone));
  } else {
    // categorical variable
    splitter_.reset(new categorical_splitter_strategy(min_num_node_obs,
						      proposedsplit_.split_class()));
  }
}


void VarSplitter::IncorporateObs(double xval, double residval, double weight) {
  splitter_->incorporate_obs(bestsplit_, proposedsplit_,
			     xval, residval, weight);
}


void VarSplitter::WrapUpCurrentVariable() {
  splitter_->wrap_up(bestsplit_, proposedsplit_);
  
  if (proposedsplit_.split_variable() == bestsplit_.split_variable()) {
    if (proposedsplit_.has_missing()) {
      bestsplit_.set_missing_def(proposedsplit_.get_missing_def());
    } else  // DEBUG: consider a weighted average with parent node?
    {
      bestsplit_.set_missing_def(NodeDef(initial_.get_weightresid(),
					 initial_.get_totalweight(),
					 0));
    }
  }
}
