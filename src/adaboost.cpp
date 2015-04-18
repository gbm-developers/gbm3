// GBM by Greg Ridgeway  Copyright (C) 2003

#include "adaboost.h"

CAdaBoost::CAdaBoost()
{
}

CAdaBoost::~CAdaBoost()
{
}


void CAdaBoost::ComputeWorkingResponse
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
  
  if(adOffset == NULL)
    {
      for(i=0; i<nTrain; i++)
        {
	  adZ[i] = -(2*adY[i]-1) * std::exp(-(2*adY[i]-1)*adF[i]);
        }
    }
  else
    {
      for(i=0; i<nTrain; i++)
        {
	  adZ[i] = -(2*adY[i]-1) * std::exp(-(2*adY[i]-1)*(adOffset[i]+adF[i]));
        }
    }
}



void CAdaBoost::InitF
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
    double dNum = 0.0;
    double dDen = 0.0;

    dInitF = 0.0;

    if(adOffset == NULL)
    {
      for(i=0; i<cLength; i++)
        {
	  if(adY[i]==1.0)
            {
	      dNum += adWeight[i];
            }
	  else
            {
                dDen += adWeight[i];
            }
        }
    }
    else
      {
        for(i=0; i<cLength; i++)
        {
	  if(adY[i]==1.0)
            {
	      dNum += adWeight[i] * std::exp(-adOffset[i]);
            }
	  else
            {
	      dDen += adWeight[i] * std::exp(adOffset[i]);
            }
        }
      }
    
    dInitF = 0.5*std::log(dNum/dDen);
}


double CAdaBoost::Deviance
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
    double dW = 0.0;

    if(adOffset == NULL)
    {
        for(i=cIdxOff; i<cLength+cIdxOff; i++)
        {
            dL += adWeight[i] * std::exp(-(2*adY[i]-1)*adF[i]);
            dW += adWeight[i];
        }
    }
    else
    {
        for(i=cIdxOff; i<cLength+cIdxOff; i++)
        {
            dL += adWeight[i] * std::exp(-(2*adY[i]-1)*(adOffset[i]+adF[i]));
            dW += adWeight[i];
       }
    }

    return dL/dW;
}


void CAdaBoost::FitBestConstant
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
  double dF = 0.0;
  unsigned long iObs = 0;
  unsigned long iNode = 0;
  vecdNum.resize(cTermNodes);
  vecdNum.assign(vecdNum.size(),0.0);
  vecdDen.resize(cTermNodes);
  vecdDen.assign(vecdDen.size(),0.0);
    

  for(iObs=0; iObs<nTrain; iObs++)
    {
      if(afInBag[iObs])
        {
	  dF = adF[iObs] + ((adOffset==NULL) ? 0.0 : adOffset[iObs]);
	  vecdNum[aiNodeAssign[iObs]] +=
	    adW[iObs]*(2*adY[iObs]-1)*std::exp(-(2*adY[iObs]-1)*dF);
	  vecdDen[aiNodeAssign[iObs]] +=
	    adW[iObs]*std::exp(-(2*adY[iObs]-1)*dF);
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
	      vecpTermNodes[iNode]->dPrediction =
		vecdNum[iNode]/vecdDen[iNode];
            }
        }
    }
}


double CAdaBoost::BagImprovement
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

            dReturnValue += adWeight[i]*
                (std::exp(-(2*adY[i]-1)*dF) -
                 std::exp(-(2*adY[i]-1)*(dF+dStepSize*adFadj[i])));
            dW += adWeight[i];
        }
    }

    return dReturnValue/dW;
}
