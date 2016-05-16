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
	cMinObsInNode = minNumObs;

	/*InitTotalWeight = 0.0;
	InitWeightResiduals = 0.0;
	InitNumObs = 0;
	fIsSplit = false;
	dLastXValue = -HUGE_VAL;
	*/
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

	double dWZ = 0.0;

	//if(fIsSplit) return;

	dWZ = dW*dZ;

	if(ISNA(dX))
	{
		proposedSplit.UpdateMissingNode(dWZ, dW);

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
		proposedSplit.UpdateLeftNode(dWZ, dW);
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

	  cFiniteMeans = proposedSplit.SetAndReturnNumGroupMeans();

	  // if only one group has a finite mean it will not consider
	  // might be all are missing so no categories enter here
	  for(i=0; (cFiniteMeans>1) && ((ULONG)i<cFiniteMeans-1); i++)
	    {


	      proposedSplit.SplitValue = (double) i;
	      proposedSplit.UpdateLeftNodeWithCat(i);
	      proposedSplit.setBestCategory();
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
}

/*void VarSplitter::SetForVariable(unsigned long iWhichVar, long cVarClasses)
{
	if(fIsSplit) return;
	//bestSplit.ResetSplitProperties(InitWeightResiduals, InitTotalWeight, InitNumObs);
	proposedSplit.ResetSplitProperties(InitWeightResiduals, InitTotalWeight, InitNumObs,
		  proposedSplit.SplitValue,	cVarClasses, iWhichVar);



}*/

void VarSplitter::ResetForNewVar
(
    unsigned long iWhichVar,
    long cCurrentVarClasses
)
{
  //if(fIsSplit) return;
  proposedSplit.ResetSplitProperties(dInitSumZ, dInitTotalW, cInitN,
  		  proposedSplit.SplitValue,	cCurrentVarClasses, iWhichVar);
  dLastXValue = -HUGE_VAL;
}

void VarSplitter::WrapUpCurrentVariable()
{
  if(proposedSplit.SplitVar == bestSplit.SplitVar)
    {
      if(proposedSplit.missing.numObs > 0)
        {
    	  bestSplit.missing.weightResid = proposedSplit.missing.weightResid;
		bestSplit.missing.totalWeight = proposedSplit.missing.totalWeight;
		bestSplit.missing.numObs = proposedSplit.missing.numObs;

        }
      else // DEBUG: consider a weighted average with parent node?
        {
		bestSplit.missing.weightResid  = dInitSumZ;
		bestSplit.missing.totalWeight = dInitTotalW;
		bestSplit.missing.numObs      = 0;
        }
    }
}

/*void VarSplitter::Reset()
{
	// Reset the splitter for new searching
	InitTotalWeight = 0.0;
	InitWeightResiduals = 0.0;
	InitNumObs = 0;

	dLastXValue = -HUGE_VAL;

	// Reset best split
	bestSplit.ResetSplitProperties(0.0, 0.0, 0);
	proposedSplit.ResetSplitProperties(0.0, 0.0, 0);

}*/

//---------------------
// Private Functions
//---------------------

/*
void VarSplitter::WrapUpSplit()
{
  if (!proposedSplit.hasMissing())
	{
		bestSplit.missing.weightResid = InitWeightResiduals;
		bestSplit.missing.totalWeight = InitTotalWeight;
		bestSplit.missing.numObs      = 0;
	}
	else
	{
	  bestSplit.missing = proposedSplit.missing;
	}
}

*/


