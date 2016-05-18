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
VarSplitter::VarSplitter(unsigned long minNumObs):bestSplit(), proposedSplit(),
adGroupSumZ(1024), adGroupW(1024), acGroupN(1024), groupMeanAndCat(1024)
{
	cMinObsInNode = minNumObs;
	fIsSplit = false;
}

VarSplitter::~VarSplitter()
{
}

void VarSplitter::IncorporateObs
(
    const double& dX,
    const double& dZ,
    const double& dW,
    const long& lMonotone
)
{
	if(fIsSplit) return;

	if(ISNA(dX))
	{
		proposedSplit.UpdateMissingNode(dW*dZ, dW);

	}
	else if(proposedSplit.SplitClass == 0)   // variable is continuous
	{
		if(dLastXValue > dX)
		{
			throw GBM::failure("Observations are not in order. gbm() was unable to build an index for the design matrix. Could be a bug in gbm or an unusual data type in data.");
		}

		// Evaluate the current split
		// the newest observation is still in the right child
		proposedSplit.SplitValue = 0.5*(dLastXValue + dX);

		if((dLastXValue != dX) &&
			proposedSplit.HasMinNumOfObs(cMinObsInNode) &&
			proposedSplit.SplitIsCorrMonotonic(lMonotone))
		{
			proposedSplit.NodeGradResiduals();
			if(proposedSplit.ImprovedResiduals > bestSplit.ImprovedResiduals)
			{
				bestSplit = proposedSplit;
			}
		}

		// now move the new observation to the left
		// if another observation arrives we will evaluate this
		proposedSplit.UpdateLeftNode(dW*dZ, dW);
		dLastXValue = dX;
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

	  if(fIsSplit) return;
	  cFiniteMeans = SetAndReturnNumGroupMeans();

	  // if only one group has a finite mean it will not consider
	  // might be all are missing so no categories enter here
	  for(i=0; (cFiniteMeans>1) && ((ULONG)i<cFiniteMeans-1); i++)
	    {


	      proposedSplit.SplitValue = (double) i;
	      UpdateLeftNodeWithCat(i);
	      proposedSplit.setBestCategory(groupMeanAndCat);
	      proposedSplit.NodeGradResiduals();

		  if(proposedSplit.HasMinNumOfObs(cMinObsInNode)
		      		  && (proposedSplit.ImprovedResiduals > bestSplit.ImprovedResiduals))
		{

		  bestSplit = proposedSplit;
	        }
	    }
}

void VarSplitter::Set(CNode& nodeToSplit)
{
	dInitSumZ = nodeToSplit.dPrediction * nodeToSplit.dTrainW;
	dInitTotalW = nodeToSplit.dTrainW;
	cInitN = nodeToSplit.cN;

	bestSplit.ResetSplitProperties(dInitSumZ, dInitTotalW, cInitN);
	fIsSplit=false;
}

void VarSplitter::ResetForNewVar
(
    unsigned long iWhichVar,
    long cCurrentVarClasses
)
{
  if(fIsSplit) return;
  proposedSplit.ResetSplitProperties(dInitSumZ, dInitTotalW, cInitN,
  		  proposedSplit.SplitValue,	cCurrentVarClasses, iWhichVar);

  std::fill(adGroupSumZ.begin(), adGroupSumZ.begin() + cCurrentVarClasses, 0);
  std::fill(adGroupW.begin(), adGroupW.begin() + cCurrentVarClasses, 0);
  std::fill(acGroupN.begin(), acGroupN.begin() + cCurrentVarClasses, 0);

  dLastXValue = -HUGE_VAL;
}

void VarSplitter::WrapUpCurrentVariable()
{
  if(proposedSplit.SplitVar == bestSplit.SplitVar)
    {
      if(proposedSplit.missing.hasObs())
	{
	  bestSplit.missing = proposedSplit.missing;
	}
      else // DEBUG: consider a weighted average with parent node?
	{
	  bestSplit.missing = NodeDef(dInitSumZ, dInitTotalW, 0);
	}
    }
}


