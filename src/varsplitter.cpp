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
VarSplitter::VarSplitter(unsigned long min_num_node_obs): bestsplit_(), proposedsplit_(),
group_sumresid_(1024), group_weight_(1024), group_num_obs_(1024), groupMeanAndCat(1024)
{
	min_num_node_obs_ = min_num_node_obs;
	issplit_ = false;
}

VarSplitter::~VarSplitter()
{
}

void VarSplitter::IncorporateObs
(
    double xval,
    double residval,
    double weight,
    long monotonicity
)
{
	if(issplit_) return;

	if(ISNA(xval))
	{
		proposedsplit_.UpdateMissingNode(weight*residval, weight);

	}
	else if(proposedsplit_.split_class_ == 0)   // variable is continuous
	{
	  if(last_xvalue_ > xval)
	    {
	      throw gbm_exception::Failure("Observations are not in order. gbm() was unable to build an index for the design matrix. Could be a bug in gbm or an unusual data type in data.");
	    }
	  
	  // Evaluate the current split
	  // the newest observation is still in the right child
	  proposedsplit_.split_value_ = 0.5*(last_xvalue_ + xval);
	  
	  if((last_xvalue_ != xval) &&
	     proposedsplit_.has_min_num_obs(min_num_node_obs_) &&
	     proposedsplit_.split_is_correct_monotonicity(monotonicity))
	    {
	      proposedsplit_.NodeGradResiduals();
	      if(proposedsplit_.improvement_ > bestsplit_.improvement_)
		{
		  bestsplit_ = proposedsplit_;
		}
	    }
	  
	  // now move the new observation to the left
	  // if another observation arrives we will evaluate this
	  proposedsplit_.UpdateLeftNode(weight*residval, weight);
	  last_xvalue_ = xval;
	}
	else // variable is categorical, evaluates later
	  {
	    IncrementCategories((unsigned long) xval, weight*residval, weight);
	  }
}


void VarSplitter::EvaluateCategoricalSplit()
{
  long i=0;
  unsigned long num_finite_means = 0;
  
  if(issplit_) return;
  num_finite_means = SetAndReturnNumGroupMeans();
  
  // if only one group has a finite mean it will not consider
  // might be all are missing so no categories enter here
  for(i=0; (num_finite_means>1) && ((unsigned long)i<num_finite_means-1); i++)
    {
      
      
      proposedsplit_.split_value_ = (double) i;
      UpdateLeftNodeWithCat(i);
      proposedsplit_.SetBestCategory(groupMeanAndCat);
      proposedsplit_.NodeGradResiduals();
      
      if(proposedsplit_.has_min_num_obs(min_num_node_obs_)
	 && (proposedsplit_.improvement_ > bestsplit_.improvement_))
	{
	  
	  bestsplit_ = proposedsplit_;
	}
    }
}

void VarSplitter::Set(CNode& node_to_split)
{
  initial_sumresiduals   = node_to_split.prediction * node_to_split.totalweight;
  initial_totalweight = node_to_split.totalweight;
  initial_numobs      = node_to_split.numobs;
  
  bestsplit_.ResetSplitProperties(initial_sumresiduals, initial_totalweight, initial_numobs);
  issplit_=false;
}

void VarSplitter::ResetForNewVar
(
 unsigned long whichvar,
 long numvar_classes
)
{
  if(issplit_) return;
  proposedsplit_.ResetSplitProperties(initial_sumresiduals,
				     initial_totalweight,
				     initial_numobs,
				     proposedsplit_.split_value_,
				     numvar_classes, whichvar);

  std::fill(group_sumresid_.begin(),
	    group_sumresid_.begin() + numvar_classes,
	    0);
  std::fill(group_weight_.begin(),
	    group_weight_.begin() + numvar_classes,
	    0);
  std::fill(group_num_obs_.begin(),
	    group_num_obs_.begin() + numvar_classes,
	    0);

  last_xvalue_ = -HUGE_VAL;
}

void VarSplitter::WrapUpCurrentVariable()
{
  if(proposedsplit_.split_var_ == bestsplit_.split_var_)
    {
      if(proposedsplit_.missing_.has_obs())
	{
	  bestsplit_.missing_ = proposedsplit_.missing_;
	}
      else // DEBUG: consider a weighted average with parent node?
	{
	  bestsplit_.missing_ = NodeDef(initial_sumresiduals, initial_totalweight, 0);
	}
    }
}


