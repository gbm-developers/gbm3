//  GBM by Greg Ridgeway  Copyright (C) 2003

#include "laplace.h"

CLaplace::CLaplace()
{
   mpLocM = NULL;
}

CLaplace::~CLaplace()
{
   if(mpLocM != NULL)
   {
      delete mpLocM;
   }
}


GBMRESULT CLaplace::ComputeWorkingResponse
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
    unsigned long i = 0;

    if(adOffset == NULL)
    {
        for(i=0; i<nTrain; i++)
        {
            adZ[i] = (adY[i] - adF[i]) > 0.0 ? 1.0 : -1.0;
        }
    }
    else
    {
        for(i=0; i<nTrain; i++)
        {
            adZ[i] = (adY[i] - adOffset[i] - adF[i]) > 0.0 ? 1.0 : -1.0;
        }
    }

    return GBM_OK;
}



GBMRESULT CLaplace::InitF
(
    double *adY,
    double *adMisc,
    double *adOffset, 
    double *adWeight,
    double &dInitF, 
    unsigned long cLength
)
{
    GBMRESULT hr = GBM_OK;

    double dOffset = 0.0;
    unsigned long ii = 0;
    int nLength = int(cLength);

    double *adArr = NULL;

    // Create a new LocationM object (for weighted medians)
    double *pTemp = NULL;
    mpLocM = new CLocationM("Other", 0, pTemp);
    if(mpLocM == NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }
    
    adArr = new double[cLength];
    if(adArr == NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }

    for (ii = 0; ii < cLength; ii++)
    {
        dOffset = (adOffset==NULL) ? 0.0 : adOffset[ii];
        adArr[ii] = adY[ii] - dOffset;
    }

    dInitF = mpLocM->Median(nLength, adArr, adWeight);

Cleanup:
    return hr;
Error:
    goto Cleanup;
}


double CLaplace::Deviance
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
            dL += adWeight[i]*fabs(adY[i]-adF[i]);
            dW += adWeight[i];
        }
    }
    else
    {
        for(i=cIdxOff; i<cLength+cIdxOff; i++)
        {
            dL += adWeight[i]*fabs(adY[i]-adOffset[i]-adF[i]);
            dW += adWeight[i];
        }
    }

    return dL/dW;
}


// DEBUG: needs weighted median
GBMRESULT CLaplace::FitBestConstant
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
    GBMRESULT hr = GBM_OK;

    unsigned long iNode = 0;
    unsigned long iObs = 0;
    unsigned long iVecd = 0;
    double dOffset;

//    vecd.resize(nTrain); // should already be this size from InitF
  
   double *adArr = new double[nTrain];
   double *adW2 = new double[nTrain];

    for(iNode=0; iNode<cTermNodes; iNode++)
    {
        if(vecpTermNodes[iNode]->cN >= cMinObsInNode)
        {
         iVecd = 0;
            for(iObs=0; iObs<nTrain; iObs++)
            {
                if(afInBag[iObs] && (aiNodeAssign[iObs] == iNode))
                {
                    dOffset = (adOffset==NULL) ? 0.0 : adOffset[iObs];
                    adArr[iVecd] = adY[iObs] - dOffset - adF[iObs];
               adW2[iVecd] = adW[iObs];
                iVecd++;
            }

            
            }

         vecpTermNodes[iNode]->dPrediction = mpLocM->Median(iVecd, adArr, adW2);

        }
    }

    return hr;
}



double CLaplace::BagImprovement
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
            
            dReturnValue += 
                adWeight[i]*(fabs(adY[i]-dF) - fabs(adY[i]-dF-dStepSize*adFadj[i]));
            dW += adWeight[i];
        }
    }

    return dReturnValue/dW;
}
