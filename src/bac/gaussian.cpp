//  GBM by Greg Ridgeway  Copyright (C) 2003

#include "gaussian.h"

CGaussian::CGaussian()
{
}

CGaussian::~CGaussian()
{
}


GBMRESULT CGaussian::ComputeWorkingResponse
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adF, 
    double *adZ, 
    double *adWeight,
    bool *afInBag,
    unsigned long nTrain,
	int cIdxOff
)
{
    GBMRESULT hr = GBM_OK;
    unsigned long i = 0;

    if((adY == NULL) || (adF == NULL) || (adZ == NULL) || (adWeight == NULL))
    {
        hr = GBM_INVALIDARG;
        goto Error;
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

Cleanup:
    return hr;
Error:
    goto Cleanup;
}



GBMRESULT CGaussian::InitF
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

    return GBM_OK;
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


GBMRESULT CGaussian::FitBestConstant
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adW,
    double *adF,
    double *adZ,
    unsigned long *aiNodeAssign,
    unsigned long nTrain,
    VEC_P_NODETERMINAL vecpTermNodes,
    unsigned long cTermNodes,
    unsigned long cMinObsInNode,
    bool *afInBag,
    double *adFadj,
	int cIdxOff
)
{
    // the tree aready stores the mean prediction
    // no refitting necessary

    return GBM_OK;
}

double CGaussian::BagImprovement
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adWeight,
    double *adF,
    double *adFadj,
    bool *afInBag,
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



