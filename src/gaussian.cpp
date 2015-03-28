//  GBM by Greg Ridgeway  Copyright (C) 2003

#include "gaussian.h"

CGaussian::CGaussian()
{
}

CGaussian::~CGaussian()
{
}


void CGaussian::ComputeWorkingResponse
(
 double *adY,
 double *adMisc,
 double *adOffset,
 double *adF,
 double *adZ,
 double *adWeight,
 int *afInBag,
 unsigned long nTrain,
 int cIdxOff
 )
{
  unsigned long i = 0;
  
  if (!(adY && adF && adZ && adWeight)) {
    throw GBM::invalid_argument();
  }
  
  if(adOffset == NULL)
    {
      for(i=0; i<nTrain; i++)
        {
	  adZ[i] = adY[i] - adF[i];
        }
    }
  else
    {
      for(i=0; i<nTrain; i++)
        {
	  adZ[i] = adY[i] - adOffset[i] - adF[i];
        }
    }
}

void CGaussian::InitF
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adWeight,
    double &dInitF,
    unsigned long cLength
)
{
    double dSum=0.0;
    double dTotalWeight = 0.0;
    unsigned long i=0;

    // compute the mean
    if(adOffset==NULL)
    {
        for(i=0; i<cLength; i++)
        {
            dSum += adWeight[i]*adY[i];
            dTotalWeight += adWeight[i];
        }
    }
    else
    {
        for(i=0; i<cLength; i++)
        {
            dSum += adWeight[i]*(adY[i] - adOffset[i]);
            dTotalWeight += adWeight[i];
        }
    }
    dInitF = dSum/dTotalWeight;
}


double CGaussian::Deviance
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adWeight,
    double *adF,
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
            dL += adWeight[i]*(adY[i]-adF[i])*(adY[i]-adF[i]);
            dW += adWeight[i];
        }
    }
    else
    {
        for(i=cIdxOff; i<cLength+cIdxOff; i++)
        {
            dL += adWeight[i]*(adY[i]-adOffset[i]-adF[i])*
                              (adY[i]-adOffset[i]-adF[i]);
            dW += adWeight[i];
       }
    }

    return dL/dW;
}


void CGaussian::FitBestConstant
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adW,
    double *adF,
    double *adZ,
    const std::vector<unsigned long>& aiNodeAssign,
    unsigned long nTrain,
    VEC_P_NODETERMINAL vecpTermNodes,
    unsigned long cTermNodes,
    unsigned long cMinObsInNode,
    int *afInBag,
    double *adFadj,
	int cIdxOff
)
{
  // the tree aready stores the mean prediction
  // no refitting necessary
}

double CGaussian::BagImprovement
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adWeight,
    double *adF,
    double *adFadj,
    int *afInBag,
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

            dReturnValue += adWeight[i]*dStepSize*adFadj[i]*
                            (2.0*(adY[i]-dF) - dStepSize*adFadj[i]);
            dW += adWeight[i];
        }
    }

    return dReturnValue/dW;
}



