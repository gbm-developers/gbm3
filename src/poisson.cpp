//  GBM by Greg Ridgeway  Copyright (C) 2003

#include "poisson.h"

CPoisson::CPoisson()
{
}

CPoisson::~CPoisson()
{
}


GBMRESULT CPoisson::ComputeWorkingResponse
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
    double dF = 0.0;

    // compute working response
    for(i=0; i < nTrain; i++)
    {
        dF = adF[i] + ((adOffset==NULL) ? 0.0 : adOffset[i]);
        adZ[i] = adY[i] - exp(dF);
    }

    return GBM_OK;
}



GBMRESULT CPoisson::InitF
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

    double dSum = 0.0;
    double dDenom = 0.0;
    unsigned long i = 0;

    if(adOffset == NULL)
    {
        for(i=0; i<cLength; i++)
        {
            dSum += adWeight[i]*adY[i];
            dDenom += adWeight[i];
        }
    }
    else
    {
        for(i=0; i<cLength; i++)
        {
            dSum += adWeight[i]*adY[i];
            dDenom += adWeight[i]*exp(adOffset[i]);
        }
    }

    dInitF = log(dSum/dDenom);

    return hr;
}


double CPoisson::Deviance
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
            dL += adWeight[i]*(adY[i]*adF[i] - exp(adF[i]));
            dW += adWeight[i];
        }
    }
    else
    {
        for(i=cIdxOff; i<cLength+cIdxOff; i++)
        {
            dL += adWeight[i]*(adY[i]*(adOffset[i]+adF[i]) - 
                               exp(adOffset[i]+adF[i]));
            dW += adWeight[i];
       }
    }

    return -2*dL/dW;
}


GBMRESULT CPoisson::FitBestConstant
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

    unsigned long iObs = 0;
    unsigned long iNode = 0;
    vecdNum.resize(cTermNodes);
    vecdNum.assign(vecdNum.size(),0.0);
    vecdDen.resize(cTermNodes);
    vecdDen.assign(vecdDen.size(),0.0);
    
    vecdMax.resize(cTermNodes);
    vecdMax.assign(vecdMax.size(),-HUGE_VAL);
    vecdMin.resize(cTermNodes);
    vecdMin.assign(vecdMin.size(),HUGE_VAL);

    if(adOffset == NULL)
    {
        for(iObs=0; iObs<nTrain; iObs++)
        {
            if(afInBag[iObs])
            {
                vecdNum[aiNodeAssign[iObs]] += adW[iObs]*adY[iObs];
                vecdDen[aiNodeAssign[iObs]] += adW[iObs]*exp(adF[iObs]);
            }
            vecdMax[aiNodeAssign[iObs]] = 
               fmax2(adF[iObs],vecdMax[aiNodeAssign[iObs]]);
            vecdMin[aiNodeAssign[iObs]] =  
               fmin2(adF[iObs],vecdMin[aiNodeAssign[iObs]]);
        }
    }
    else
    {
        for(iObs=0; iObs<nTrain; iObs++)
        {
            if(afInBag[iObs])
            {
                vecdNum[aiNodeAssign[iObs]] += adW[iObs]*adY[iObs];
                vecdDen[aiNodeAssign[iObs]] += 
                    adW[iObs]*exp(adOffset[iObs]+adF[iObs]);
            }
        }        
    }
    for(iNode=0; iNode<cTermNodes; iNode++)
    {
        if(vecpTermNodes[iNode]!=NULL)
        {
            if(vecdNum[iNode] == 0.0)
            {
                // DEBUG: if vecdNum==0 then prediction = -Inf
                // Not sure what else to do except plug in an arbitrary
                //   negative number, -1? -10? Let's use -1, then make
                //   sure |adF| < 19 always.
                vecpTermNodes[iNode]->dPrediction = -19.0;
            }
            else if(vecdDen[iNode] == 0.0)
            {
                vecpTermNodes[iNode]->dPrediction = 0.0;
            }            
            else
            {
                vecpTermNodes[iNode]->dPrediction = 
                    log(vecdNum[iNode]/vecdDen[iNode]);
            }
            vecpTermNodes[iNode]->dPrediction = 
               fmin2(vecpTermNodes[iNode]->dPrediction,
                     19-vecdMax[iNode]);
            vecpTermNodes[iNode]->dPrediction = 
               fmax2(vecpTermNodes[iNode]->dPrediction,
                     -19-vecdMin[iNode]);
        }
    }

    return hr;
}


double CPoisson::BagImprovement
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

            dReturnValue += adWeight[i]*
                            (adY[i]*dStepSize*adFadj[i] -
                             exp(dF+dStepSize*adFadj[i]) +
                             exp(dF));
            dW += adWeight[i];
        }
    }

    return dReturnValue/dW;
}


