// GBM by Greg Ridgeway  Copyright (C) 2003

#include "bernoulli.h"

CBernoulli::CBernoulli()
{
  // Used to issue warnings to user that at least one terminal node capped
  fCappedPred = false;
}

CBernoulli::~CBernoulli()
{
}


void CBernoulli::ComputeWorkingResponse
(
    const double *adY,
    const double *adMisc,
    const double *adOffset,
    const double *adF,
    double *adZ,
    const double *adWeight,
    const bag& afInBag,
    unsigned long nTrain,
    int cIdxOff
)
{
  unsigned long i = 0;
  double dProb = 0.0;
  double dF = 0.0;

  for(i=0; i<nTrain; i++)
  {
    dF = adF[i] + ((adOffset==NULL) ? 0.0 : adOffset[i]);
    dProb = 1.0/(1.0+std::exp(-dF));

    adZ[i] = adY[i] - dProb;
#ifdef NOISY_DEBUG
//  Rprintf("dF=%f, dProb=%f, adZ=%f, adY=%f\n", dF, dProb, adZ[i], adY[i]);
    if(dProb<  0.0001) Rprintf("Small prob(i=%d)=%f Z=%f\n",i,dProb,adZ[i]);
    if(dProb>1-0.0001) Rprintf("Large prob(i=%d)=%f Z=%f\n",i,dProb,adZ[i]);
#endif
  }
}


void CBernoulli::InitF
(
    const double *adY,
    const double *adMisc,
    const double *adOffset,
    const double *adWeight,
    double &dInitF,
    unsigned long cLength
)
{
    unsigned long i=0;
    double dTemp=0.0;

    if(adOffset==NULL)
    {
        double dSum=0.0;
        for(i=0; i<cLength; i++)
        {
            dSum += adWeight[i]*adY[i];
            dTemp += adWeight[i];
        }
        dInitF = std::log(dSum/(dTemp-dSum));
    }
    else
    {
        // Newton method for solving for F
        // should take about 3-6 iterations.
        double dNum=0.0;         // numerator
        double dDen=0.0;         // denominator
        double dNewtonStep=1.0;  // change
        dInitF = 0.0;
        while(dNewtonStep > 0.0001)
        {
            dNum=0.0;
            dDen=0.0;
            for(i=0; i<cLength; i++)
            {
                dTemp = 1.0/(1.0+std::exp(-(adOffset[i] + dInitF)));
                dNum += adWeight[i]*(adY[i]-dTemp);
                dDen += adWeight[i]*dTemp*(1.0-dTemp);
            }
            dNewtonStep = dNum/dDen;
            dInitF += dNewtonStep;
        }
    }
}



double CBernoulli::Deviance
(
    const double *adY,
    const double *adMisc,
    const double *adOffset,
    const double *adWeight,
    const double *adF,
    unsigned long cLength,
    int cIdxOff
)
{
   unsigned long i=0;
   double dL = 0.0;
   double dF = 0.0;
   double dW = 0.0;

   if(adOffset==NULL)
   {
      for(i=cIdxOff; i<cLength+cIdxOff; i++)
      {
         dL += adWeight[i]*(adY[i]*adF[i] - std::log(1.0+std::exp(adF[i])));
         dW += adWeight[i];
      }
   }
   else
   {
      for(i=cIdxOff; i<cLength+cIdxOff; i++)
      {
         dF = adF[i] + adOffset[i];
         dL += adWeight[i]*(adY[i]*dF - std::log(1.0+std::exp(dF)));
         dW += adWeight[i];
      }
   }

   return -2*dL/dW;
}


void CBernoulli::FitBestConstant
(
  const double *adY,
  const double *adMisc,
  const double *adOffset,
  const double *adW,
  const double *adF,
  double *adZ,
  const std::vector<unsigned long>& aiNodeAssign,
  unsigned long nTrain,
  VEC_P_NODETERMINAL vecpTermNodes,
  unsigned long cTermNodes,
  unsigned long cMinObsInNode,
  const bag& afInBag,
  const double *adFadj,
  int cIdxOff
)
{
  unsigned long iObs = 0;
  unsigned long iNode = 0;
  double dTemp = 0.0;
  
  vecdNum.resize(cTermNodes);
  vecdNum.assign(vecdNum.size(),0.0);
  vecdDen.resize(cTermNodes);
  vecdDen.assign(vecdDen.size(),0.0);

  for(iObs=0; iObs<nTrain; iObs++)
  {
    if(afInBag[iObs])
    {
      vecdNum[aiNodeAssign[iObs]] += adW[iObs]*adZ[iObs];
      vecdDen[aiNodeAssign[iObs]] +=
          adW[iObs]*(adY[iObs]-adZ[iObs])*(1-adY[iObs]+adZ[iObs]);
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
    if(vecpTermNodes[iNode]!=NULL)
    {
      if(vecdDen[iNode] == 0)
      {
          vecpTermNodes[iNode]->dPrediction = 0.0;
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
        vecpTermNodes[iNode]->dPrediction = dTemp;              
      }
    }
  }
}


double CBernoulli::BagImprovement
(
    const double *adY,
    const double *adMisc,
    const double *adOffset,
    const double *adWeight,
    const double *adF,
    const double *adFadj,
    const bag& afInBag,
    double dStepSize,
    unsigned long nTrain
)
{
    double dReturnValue = 0.0;
    double dF = 0.0;
    double dW = 0.0;
    unsigned long i = 0;

    for(i=0; i<nTrain; i++)
    {
        if(!afInBag[i])
        {
            dF = adF[i] + ((adOffset==NULL) ? 0.0 : adOffset[i]);

            if(adY[i]==1.0)
            {
                dReturnValue += adWeight[i]*dStepSize*adFadj[i];
            }
            dReturnValue += adWeight[i]*
                            (std::log(1.0+std::exp(dF)) -
                             std::log(1.0+std::exp(dF+dStepSize*adFadj[i])));
            dW += adWeight[i];
        }
    }

    return dReturnValue/dW;
}
