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
VarSplitter::VarSplitter(unsigned long minNumObs):bestSplit(), proposedSplit()
{

	InitTotalWeight = 0.0;
	InitWeightResiduals = 0.0;
	InitNumObs = 0;

	dLastXValue = -HUGE_VAL;
	minObsInNode = minNumObs;
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
	if(ISNA(dX))
	{
		proposedSplit.UpdateMissingNode(dW*dZ, dW);
	}
	else if(proposedSplit.SplitClass == 0)
	{
		if(dLastXValue > dX)
		{
			throw GBM::failure("Observations are not in order. gbm() was unable to build an index for the design matrix. Could be a bug in gbm or an unusual data type in data.");
		}

		// Evaluate the current split
		// the newest observation is still in the right child
		proposedSplit.SplitValue = 0.5*(dLastXValue + dX);

		if((dLastXValue != dX) &&
			proposedSplit.HasMinNumOfObs(minObsInNode) &&
			proposedSplit.SplitIsCorrMonotonic(lMonotone))
		{
			proposedSplit.NodeGradResiduals();

			if(proposedSplit.HasMinNumOfObs(minObsInNode) &&
					(proposedSplit.ImprovedResiduals > bestSplit.ImprovedResiduals))
			{
				bestSplit = proposedSplit;
				WrapUpSplit();

			}

		}

		// now move the new observation to the left
		// if another observation arrives we will evaluate this
		proposedSplit.UpdateLeftNode(dW*dZ, dW);
		dLastXValue = dX;
	}
	else // variable is categorical, evaluates later
	{
		proposedSplit.IncrementCategories((unsigned long) dX, dW*dZ, dW);
	}
}


void VarSplitter::EvaluateCategoricalSplit()
{
  long i=0;
  unsigned long cFiniteMeans = 0;

  if(proposedSplit.SplitClass == 0)
	{
	  throw GBM::invalid_argument();
	}

  cFiniteMeans = proposedSplit.SetAndReturnNumGroupMeans();

  // if only one group has a finite mean it will not consider
  // might be all are missing so no categories enter here
  for(i=0; (cFiniteMeans>1) && ((ULONG)i<cFiniteMeans-1); i++)
    {
      proposedSplit.SplitValue = (double)i;
      proposedSplit.UpdateLeftNodeWithCat(i);
      proposedSplit.NodeGradResiduals();
      proposedSplit.setBestCategory();

      if(proposedSplit.HasMinNumOfObs(minObsInNode)
    		  && (proposedSplit.ImprovedResiduals > bestSplit.ImprovedResiduals))
      {
    	  bestSplit = proposedSplit;
		  WrapUpSplit();

      }

    }
}

void VarSplitter::SetForNode(CNode& nodeToSplit)
{
	InitWeightResiduals = nodeToSplit.dPrediction * nodeToSplit.dTrainW;
	InitTotalWeight = nodeToSplit.dTrainW;
	InitNumObs = nodeToSplit.cN;
}

void VarSplitter::SetForVariable(unsigned long iWhichVar, long cVarClasses)
{

	bestSplit.ResetSplitProperties(InitWeightResiduals, InitTotalWeight, InitNumObs);
	proposedSplit.ResetSplitProperties(InitWeightResiduals, InitTotalWeight, InitNumObs,
		  proposedSplit.SplitValue,	cVarClasses, iWhichVar);

	dLastXValue = -HUGE_VAL;

}

void VarSplitter::Reset()
{
	// Reset the splitter for new searching
	InitTotalWeight = 0.0;
	InitWeightResiduals = 0.0;
	InitNumObs = 0;

	dLastXValue = -HUGE_VAL;

	// Reset best split
	bestSplit.ResetSplitProperties(0.0, 0.0, 0);
	proposedSplit.ResetSplitProperties(0.0, 0.0, 0);

}

//---------------------
// Private Functions
//---------------------
void VarSplitter::WrapUpSplit()
{
	if(proposedSplit.MissingNumObs <= 0)
	{
		bestSplit.MissingWeightResiduals   = InitWeightResiduals;
		bestSplit.MissingTotalWeight = InitTotalWeight;
		bestSplit.MissingNumObs      = 0;
	}
	else
	{
		bestSplit.MissingWeightResiduals = proposedSplit.MissingWeightResiduals;
		bestSplit.MissingTotalWeight = proposedSplit.MissingTotalWeight;
		bestSplit.MissingNumObs = proposedSplit.MissingNumObs;

	}
}


