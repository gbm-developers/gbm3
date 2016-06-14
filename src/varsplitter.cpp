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

//---------------------
// Public Functions
//---------------------
VarSplitter::VarSplitter()
    : bestsplit_(),
      proposedsplit_(),
      group_sumresid_(1024),
      group_weight_(1024),
      group_num_obs_(1024),
      groupMeanAndCat(1024) {
  min_num_node_obs_ = 0;
  last_xvalue_ = -HUGE_VAL;
  issplit_ = false;
}

VarSplitter::VarSplitter(CNode& nodetosplit,
		unsigned long min_num_node_obs,
		unsigned long whichvar, unsigned long numvar_classes)
    : initial_sumresiduals_(nodetosplit.get_prediction() * nodetosplit.get_totalweight()),
      initial_totalweight_(nodetosplit.get_totalweight()),
      initial_numobs_(nodetosplit.get_numobs()),
      bestsplit_(initial_sumresiduals_, initial_totalweight_, initial_numobs_),
      proposedsplit_(initial_sumresiduals_, initial_totalweight_, initial_numobs_,
    		  numvar_classes, whichvar),
      group_sumresid_(1024),
      group_weight_(1024),
      group_num_obs_(1024),
      groupMeanAndCat(1024) {

  min_num_node_obs_ = min_num_node_obs;
  last_xvalue_ = -HUGE_VAL;
  issplit_ = nodetosplit.is_split_determined();
}

VarSplitter::~VarSplitter() {}

void VarSplitter::IncorporateObs(double xval, double residval, double weight,
                                 long monotonicity) {
  if (issplit_) return;
  if (ISNA(xval)) {
    proposedsplit_.UpdateMissingNode(weight * residval, weight);

  } else if (proposedsplit_.split_class() == 0)  // variable is continuous
  {
    if (last_xvalue_ > xval) {
      throw gbm_exception::Failure(
          "Observations are not in order. gbm() was unable to build an index "
          "for the design matrix. Could be a bug in gbm or an unusual data "
          "type in data.");
    }

    // Evaluate the current split
    // the newest observation is still in the right child
    proposedsplit_.set_split_value(0.5 * (last_xvalue_ + xval));

    if ((last_xvalue_ != xval) &&
        proposedsplit_.has_min_num_obs(min_num_node_obs_) &&
        proposedsplit_.split_is_correct_monotonicity(monotonicity)) {
      proposedsplit_.NodeGradResiduals();
      if (proposedsplit_.get_improvement() > bestsplit_.get_improvement()) {
        bestsplit_ = proposedsplit_;
      }
    }

    // now move the new observation to the left
    // if another observation arrives we will evaluate this
    proposedsplit_.UpdateLeftNode(weight * residval, weight);
    last_xvalue_ = xval;
  } else  // variable is categorical, evaluates later
  {
    IncrementCategories((unsigned long)xval, weight * residval, weight);
  }
}

void VarSplitter::EvaluateCategoricalSplit() {
  long i = 0;
  unsigned long num_finite_means = 0;

  if (issplit_) return;
  num_finite_means = SetAndReturnNumGroupMeans();

  // if only one group has a finite mean it will not consider
  // might be all are missing so no categories enter here
  for (i = 0;
       (num_finite_means > 1) && ((unsigned long)i < num_finite_means - 1);
       i++) {
    proposedsplit_.set_split_value((double)i);
    UpdateLeftNodeWithCat(i);
    proposedsplit_.SetBestCategory(groupMeanAndCat);
    proposedsplit_.NodeGradResiduals();

    if (proposedsplit_.has_min_num_obs(min_num_node_obs_) &&
        (proposedsplit_.get_improvement() > bestsplit_.get_improvement())) {
      bestsplit_ = proposedsplit_;
    }
  }
}

void VarSplitter::Set(CNode& nodetosplit) {
  initial_sumresiduals_ =
		  nodetosplit.get_prediction() * nodetosplit.get_totalweight();
  initial_totalweight_ = nodetosplit.get_totalweight();
  initial_numobs_ = nodetosplit.get_numobs();

  bestsplit_.ResetSplitProperties(initial_sumresiduals_, initial_totalweight_,
                                    initial_numobs_);
}

void VarSplitter::WrapUpCurrentVariable() {
  if (proposedsplit_.split_variable() == bestsplit_.split_variable()) {
    if (proposedsplit_.has_missing()) {
      bestsplit_.set_missing_def(proposedsplit_.get_missing_def());
    } else  // DEBUG: consider a weighted average with parent node?
    {
      bestsplit_.set_missing_def(
          NodeDef(initial_sumresiduals_, initial_totalweight_, 0));
    }
  }
}
