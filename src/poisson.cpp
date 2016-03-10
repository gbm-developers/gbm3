//-----------------------------------
//
// File: poisson.cpp
//
// Description: poisson distribution for GBM.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "poisson.h"

//----------------------------------------
// Function Members - Private
//----------------------------------------
CPoisson::CPoisson(SEXP radMisc): CDistribution(radMisc)
{
}


//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CPoisson::Create(SEXP radMisc,
										const char* szIRMeasure,
										int& cGroups, int& cTrain)
{
	return new CPoisson(radMisc);
}


CPoisson::~CPoisson()
{
}


void CPoisson::ComputeWorkingResponse
(
	const CDataset* pData,
    const double *adF,
    double *adZ,
    const bag& afInBag,
    unsigned long nTrain
)
{
    unsigned long i = 0;
    double dF = 0.0;

    // compute working response
    for(i=0; i < nTrain; i++)
    {
        dF = adF[i] + ((pData->offset_ptr(false)==NULL) ? 0.0 : pData->offset_ptr(false)[i]);
        adZ[i] = pData->y_ptr()[i] - std::exp(dF);
    }
}



void CPoisson::InitF
(
	const CDataset* pData,
    double &dInitF,
    unsigned long cLength
)
{
    double dSum = 0.0;
    double dDenom = 0.0;
    unsigned long i = 0;

    if(pData->offset_ptr(false) == NULL)
    {
        for(i=0; i<cLength; i++)
        {
            dSum += pData->weight_ptr()[i]*pData->y_ptr()[i];
            dDenom += pData->weight_ptr()[i];
        }
    }
    else
    {
        for(i=0; i<cLength; i++)
        {
            dSum += pData->weight_ptr()[i]*pData->y_ptr()[i];
            dDenom += pData->weight_ptr()[i]*std::exp(pData->offset_ptr(false)[i]);
        }
    }

    dInitF = std::log(dSum/dDenom);
}


double CPoisson::Deviance
(
	const CDataset* pData,
    const double *adF,
    unsigned long cLength,
    bool isValidationSet
)
{
    unsigned long i=0;
    double dL = 0.0;
    double dW = 0.0;

    // Switch to validation set if necessary
    if(isValidationSet)
    {
 	   pData->shift_to_validation();
    }

    if(pData->offset_ptr(false) == NULL)
    {
        for(i=0; i<cLength; i++)
        {
            dL += pData->weight_ptr()[i]*(pData->y_ptr()[i]*adF[i] - std::exp(adF[i]));
            dW += pData->weight_ptr()[i];
        }
    }
    else
    {
        for(i=0; i<cLength; i++)
        {
            dL += pData->weight_ptr()[i]*(pData->y_ptr()[i]*(pData->offset_ptr(false)[i]+adF[i]) -
                               std::exp(pData->offset_ptr(false)[i]+adF[i]));
            dW += pData->weight_ptr()[i];
       }
    }

    // Switch back to training set if necessary
    if(isValidationSet)
    {
 	   pData->shift_to_train();
    }

    return -2*dL/dW;
}


void CPoisson::FitBestConstant
(
	const CDataset* pData,
    const double *adF,
    double *adZ,
    const std::vector<unsigned long>& aiNodeAssign,
    unsigned long nTrain,
    VEC_P_NODETERMINAL vecpTermNodes,
    unsigned long cTermNodes,
    unsigned long cMinObsInNode,
    const bag& afInBag,
    const double *adFadj
)
{
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

    if(pData->offset_ptr(false) == NULL)
    {
        for(iObs=0; iObs<nTrain; iObs++)
        {
            if(afInBag[iObs])
            {
                vecdNum[aiNodeAssign[iObs]] += pData->weight_ptr()[iObs]*pData->y_ptr()[iObs];
                vecdDen[aiNodeAssign[iObs]] += pData->weight_ptr()[iObs]*std::exp(adF[iObs]);
            }
            vecdMax[aiNodeAssign[iObs]] =
               R::fmax2(adF[iObs],vecdMax[aiNodeAssign[iObs]]);
            vecdMin[aiNodeAssign[iObs]] =
               R::fmin2(adF[iObs],vecdMin[aiNodeAssign[iObs]]);
        }
    }
    else
    {
        for(iObs=0; iObs<nTrain; iObs++)
        {
            if(afInBag[iObs])
            {
                vecdNum[aiNodeAssign[iObs]] += pData->weight_ptr()[iObs]*pData->y_ptr()[iObs];
                vecdDen[aiNodeAssign[iObs]] +=
                    pData->weight_ptr()[iObs]*std::exp(pData->offset_ptr(false)[iObs]+adF[iObs]);
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
                    std::log(vecdNum[iNode]/vecdDen[iNode]);
            }
            vecpTermNodes[iNode]->dPrediction =
               R::fmin2(vecpTermNodes[iNode]->dPrediction,
                     19-vecdMax[iNode]);
            vecpTermNodes[iNode]->dPrediction =
               R::fmax2(vecpTermNodes[iNode]->dPrediction,
                     -19-vecdMin[iNode]);
        }
    }
}


double CPoisson::BagImprovement
(
	const CDataset* pData,
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
            dF = adF[i] + ((pData->offset_ptr(false)==NULL) ? 0.0 : pData->offset_ptr(false)[i]);

            dReturnValue += pData->weight_ptr()[i]*
                            (pData->y_ptr()[i]*dStepSize*adFadj[i] -
                             std::exp(dF+dStepSize*adFadj[i]) +
                             std::exp(dF));
            dW += pData->weight_ptr()[i];
        }
    }

    return dReturnValue/dW;
}


