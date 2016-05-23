//-----------------------------------
//
// File: bernoulli.cpp
//
// Description: bernoulli distribution.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------

#include "bernoulli.h"
#include <memory>

//----------------------------------------
// Function Members - Private
//----------------------------------------
CBernoulli::CBernoulli()
{
  // Used to issue warnings to user that at least one terminal node capped
  fCappedPred = false;
}

//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CBernoulli::Create(DataDistParams& distParams)
{
	return new CBernoulli();
}

CBernoulli::~CBernoulli()
{
}

void CBernoulli::ComputeWorkingResponse
(
	const CDataset& data,
    const double *adF,
    double *adZ
)
{
  double dProb = 0.0;
  double dF = 0.0;

  for(unsigned long i=0; i<data.get_trainSize(); i++)
  {
    dF = adF[i] +  data.offset_ptr()[i];
    dProb = 1.0/(1.0+std::exp(-dF));

    adZ[i] = data.y_ptr()[i] - dProb;
#ifdef NOISY_DEBUG
//  Rprintf("dF=%f, dProb=%f, adZ=%f, data.y_ptr()=%f\n", dF, dProb, adZ[i], data.y_ptr()[i]);
    if(dProb<  0.0001) Rprintf("Small prob(i=%d)=%f Z=%f\n",i,dProb,adZ[i]);
    if(dProb>1-0.0001) Rprintf("Large prob(i=%d)=%f Z=%f\n",i,dProb,adZ[i]);
#endif
  }
}


double CBernoulli::InitF
(
	const CDataset& data
)
{
    // Newton method for solving for F
    // should take about 3-6 iterations.

    double dInitF = 0.0;
    double dNewtonStep=1.0; // set to 1 initially

    for (int itCount=0;
	 (itCount < 6) && (std::abs(dNewtonStep) > 0.001);
	 ++itCount) {

      double dNum=0.0;
      double dDen=0.0;

      for(unsigned long i=0; i<data.get_trainSize(); i++)
	{
	  const double dTemp = 1.0/(1.0+std::exp(-(data.offset_ptr()[i] + dInitF)));
	  dNum += data.weight_ptr()[i]*(data.y_ptr()[i]-dTemp);
	  dDen += data.weight_ptr()[i]*dTemp*(1.0-dTemp);
	}
      dNewtonStep = dNum/dDen;
      dInitF += dNewtonStep;
    }


    return dInitF;
}

double CBernoulli::Deviance
(
	const CDataset& data,
    const double *adF,
    bool isValidationSet
)
{
   unsigned long i=0;
   double dL = 0.0;
   double dF = 0.0;
   double dW = 0.0;

   // Switch to validation set if necessary
   unsigned long cLength = data.get_trainSize();
   if(isValidationSet)
   {
	   data.shift_to_validation();
	   cLength = data.GetValidSize();
   }


	for(i=0; i!=cLength; i++)
	{
	 dF = adF[i] + data.offset_ptr()[i];
	 dL += data.weight_ptr()[i]*(data.y_ptr()[i]*dF - std::log(1.0+std::exp(dF)));
	 dW += data.weight_ptr()[i];
	}

   // Switch back to trainig set if necessary
   if(isValidationSet)
   {
	   data.shift_to_train();
   }

   //TODO: Check if weights are all zero for validation set
	if((dW == 0.0) && (dL == 0.0))
	{
		return nan("");
	}
	else if(dW == 0.0)
	{
		return copysign(HUGE_VAL, -dL);
	}

   return -2*dL/dW;
}


void CBernoulli::FitBestConstant
(
  const CDataset& data,
  const double *adF,
  unsigned long cTermNodes,
  double* adZ,
  CTreeComps& treeComps
)
{
  unsigned long iObs = 0;
  unsigned long iNode = 0;
  double dTemp = 0.0;
  vector<double> vecdNum(cTermNodes, 0.0);
  vector<double> vecdDen(cTermNodes, 0.0);

  for(iObs=0; iObs<data.get_trainSize(); iObs++)
  {
    if(data.GetBagElem(iObs))
    {
      vecdNum[treeComps.GetNodeAssign()[iObs]] += data.weight_ptr()[iObs]*adZ[iObs];
      vecdDen[treeComps.GetNodeAssign()[iObs]] +=
          data.weight_ptr()[iObs]*(data.y_ptr()[iObs]-adZ[iObs])*(1-data.y_ptr()[iObs]+adZ[iObs]);
#ifdef NOISY_DEBUG
/*
      Rprintf("iNode=%d, dNum(%d)=%f, dDen(%d)=%f\n",
              aiNodeAssign[iObs],
              iObs,vecdNum[aiNodeAssign[iObs]],
              iObs,vecdDen[aiNodeAssign[iObs]]);
*/
#endif
    }
  }

  for(iNode=0; iNode<cTermNodes; iNode++)
  {
    if(treeComps.GetTermNodes()[iNode]!=NULL)
    {
      if(vecdDen[iNode] == 0)
      {
          treeComps.GetTermNodes()[iNode]->dPrediction = 0.0;
      }
      else
      {
        dTemp = vecdNum[iNode]/vecdDen[iNode];
        // avoid large changes in predictions on log odds scale
        if(std::abs(dTemp) > 1.0)
        {
          if(!fCappedPred)
          {
            // set fCappedPred=true so that warning only issued once
            fCappedPred = true;  
            Rcpp::warning("Some terminal node predictions were excessively large for Bernoulli and have been capped at 1.0. Likely due to a feature that separates the 0/1 outcomes. Consider reducing shrinkage parameter.");
          }
          if(dTemp>1.0) dTemp = 1.0;
          else if(dTemp<-1.0) dTemp = -1.0;
        }
        treeComps.GetTermNodes()[iNode]->dPrediction = dTemp;
      }
    }
  }
}


double CBernoulli::BagImprovement
(
    const CDataset& data,
    const double *adF,
    const double shrinkage,
    const double* adFadj
)
{
    double dReturnValue = 0.0;
    double dF = 0.0;
    double dW = 0.0;
    unsigned long i = 0;

    for(i=0; i<data.get_trainSize(); i++)
    {
        if(!data.GetBagElem(i))
        {
            dF = adF[i] +  data.offset_ptr()[i];

            if(data.y_ptr()[i]==1.0)
            {
                dReturnValue += data.weight_ptr()[i]*shrinkage*adFadj[i];
            }
            dReturnValue += data.weight_ptr()[i]*
                            (std::log(1.0+std::exp(dF)) -
                             std::log(1.0+std::exp(dF+shrinkage*adFadj[i])));
            dW += data.weight_ptr()[i];
        }
    }

    return dReturnValue/dW;
}
