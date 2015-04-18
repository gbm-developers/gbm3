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
    const double *adY,
    const double *adMisc,
    const double *adOffset,
    const double *adWeight,
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
  // the tree aready stores the mean prediction
  // no refitting necessary
}

double CGaussian::BagImprovement
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

            dReturnValue += adWeight[i]*dStepSize*adFadj[i]*
                            (2.0*(adY[i]-dF) - dStepSize*adFadj[i]);
            dW += adWeight[i];
        }
    }

    return dReturnValue/dW;
}



