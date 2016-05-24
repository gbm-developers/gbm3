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
VarSplitter::VarSplitter(unsigned long minNumObs): bestsplit_(), proposedsplit_(),
group_sumresid_(1024), group_weight_(1024), group_num_obs_(1024), groupMeanAndCat(1024)
{
	min_num_node_obs_ = minNumObs;
	issplit_ = false;
}

VarSplitter::~VarSplitter()
{
}

void VarSplitter::IncorporateObs
(
    double dX,
    double dZ,
    double dW,
    long lMonotone
)
{
	if(issplit_) return;

	if(ISNA(dX))
	{
		proposedsplit_.UpdateMissingNode(dW*dZ, dW);

	}
	else if(proposedsplit_.split_class_ == 0)   // variable is continuous
	{
	  if(last_xvalue_ > dX)
	    {
	      throw GBM::Failure("Observations are not in order. gbm() was unable to build an index for the design matrix. Could be a bug in gbm or an unusual data type in data.");
	    }
	  
	  // Evaluate the current split
	  // the newest observation is still in the right child
	  proposedsplit_.split_value_ = 0.5*(last_xvalue_ + dX);
	  
	  if((last_xvalue_ != dX) &&
	     proposedsplit_.has_min_num_obs(min_num_node_obs_) &&
	     proposedsplit_.split_is_correct_monotonicity(lMonotone))
	    {
	      proposedsplit_.NodeGradResiduals();
	      if(proposedsplit_.improvement_ > bestsplit_.improvement_)
		{
		  bestsplit_ = proposedsplit_;
		}
	    }
	  
	  // now move the new observation to the left
	  // if another observation arrives we will evaluate this
	  proposedsplit_.UpdateLeftNode(dW*dZ, dW);
	  last_xvalue_ = dX;
	}
	else // variable is categorical, evaluates later
	  {
	    IncrementCategories((unsigned long) dX, dW*dZ, dW);
	  }
}


void VarSplitter::EvaluateCategoricalSplit()
{
  long i=0;
  unsigned long cFiniteMeans = 0;
  
  if(issplit_) return;
  cFiniteMeans = SetAndReturnNumGroupMeans();
  
  // if only one group has a finite mean it will not consider
  // might be all are missing so no categories enter here
  for(i=0; (cFiniteMeans>1) && ((ULONG)i<cFiniteMeans-1); i++)
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

void VarSplitter::Set(CNode& nodeToSplit)
{
  initial_sumresiduals   = nodeToSplit.prediction * nodeToSplit.totalweight;
  initial_totalweight = nodeToSplit.totalweight;
  initial_numobs      = nodeToSplit.numobs;
  
  bestsplit_.ResetSplitProperties(initial_sumresiduals, initial_totalweight, initial_numobs);
  issplit_=false;
}

void VarSplitter::ResetForNewVar
(
 unsigned long iWhichVar,
 long cCurrentVarClasses
)
{
  if(issplit_) return;
  proposedsplit_.ResetSplitProperties(initial_sumresiduals,
				     initial_totalweight,
				     initial_numobs,
				     proposedsplit_.split_value_,
				     cCurrentVarClasses, iWhichVar);

  std::fill(group_sumresid_.begin(),
	    group_sumresid_.begin() + cCurrentVarClasses,
	    0);
  std::fill(group_weight_.begin(),
	    group_weight_.begin() + cCurrentVarClasses,
	    0);
  std::fill(group_num_obs_.begin(),
	    group_num_obs_.begin() + cCurrentVarClasses,
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


