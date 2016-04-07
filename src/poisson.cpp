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
								int& cTrain)
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
    double *adZ
)
{
    unsigned long i = 0;
    double dF = 0.0;

    // compute working response
    for(i=0; i < pData->get_trainSize(); i++)
    {
        dF = adF[i] +  pData->offset_ptr(false)[i];
        adZ[i] = pData->y_ptr()[i] - std::exp(dF);
    }
}



double CPoisson::InitF
(
	const CDataset* pData
)
{
    double dSum = 0.0;
    double dDenom = 0.0;
    unsigned long i = 0;


	for(i=0; i<pData->get_trainSize(); i++)
	{
		dSum += pData->weight_ptr()[i]*pData->y_ptr()[i];
		dDenom += pData->weight_ptr()[i]*std::exp(pData->offset_ptr(false)[i]);
	}

    return std::log(dSum/dDenom);
}


double CPoisson::Deviance
(
	const CDataset* pData,
    const double *adF,
    bool isValidationSet
)
{
    unsigned long i=0;
    double dL = 0.0;
    double dW = 0.0;

    // Switch to validation set if necessary
    long cLength = pData->get_trainSize();
    if(isValidationSet)
    {
 	   pData->shift_to_validation();
 	   cLength = pData->GetValidSize();
    }


	for(i=0; i<cLength; i++)
	{
		dL += pData->weight_ptr()[i]*(pData->y_ptr()[i]*(pData->offset_ptr(false)[i]+adF[i]) -
						   std::exp(pData->offset_ptr(false)[i]+adF[i]));
		dW += pData->weight_ptr()[i];
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
    unsigned long cTermNodes,
    double* adZ,
    CTreeComps* pTreeComps
)
{
    unsigned long iObs = 0;
    unsigned long iNode = 0;
    vector<double> vecdNum(cTermNodes, 0.0);
    vector<double> vecdDen(cTermNodes, 0.0);
    vector<double> vecdMax(cTermNodes, -HUGE_VAL);
    vector<double> vecdMin(cTermNodes, HUGE_VAL);

	for(iObs=0; iObs<pData->get_trainSize(); iObs++)
	{
		if(pData->GetBagElem(iObs))
		{
			vecdNum[pTreeComps->GetNodeAssign()[iObs]] += pData->weight_ptr()[iObs]*pData->y_ptr()[iObs];
			vecdDen[pTreeComps->GetNodeAssign()[iObs]] +=
				pData->weight_ptr()[iObs]*std::exp(pData->offset_ptr(false)[iObs]+adF[iObs]);
		}
	}

    for(iNode=0; iNode<cTermNodes; iNode++)
    {
        if(pTreeComps->GetTermNodes()[iNode]!=NULL)
        {
            if(vecdNum[iNode] == 0.0)
            {
                // DEBUG: if vecdNum==0 then prediction = -Inf
                // Not sure what else to do except plug in an arbitrary
                //   negative number, -1? -10? Let's use -1, then make
                //   sure |adF| < 19 always.
            	pTreeComps->GetTermNodes()[iNode]->dPrediction = -19.0;
            }
            else if(vecdDen[iNode] == 0.0)
            {
            	pTreeComps->GetTermNodes()[iNode]->dPrediction = 0.0;
            }
            else
            {
            	pTreeComps->GetTermNodes()[iNode]->dPrediction =
                    std::log(vecdNum[iNode]/vecdDen[iNode]);
            }
            pTreeComps->GetTermNodes()[iNode]->dPrediction =
               R::fmin2(pTreeComps->GetTermNodes()[iNode]->dPrediction,
                     19-vecdMax[iNode]);
            pTreeComps->GetTermNodes()[iNode]->dPrediction =
               R::fmax2(pTreeComps->GetTermNodes()[iNode]->dPrediction,
                     -19-vecdMin[iNode]);
        }
    }
}


double CPoisson::BagImprovement
(
	const CDataset& data,
    const double *adF,
    const bag& afInBag,
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
            dF = adF[i] + data.offset_ptr(false)[i];

            dReturnValue += data.weight_ptr()[i]*
                            (data.y_ptr()[i]*shrinkage*adFadj[i] -
                             std::exp(dF+shrinkage*adFadj[i]) +
                             std::exp(dF));
            dW += data.weight_ptr()[i];
        }
    }

    return dReturnValue/dW;
}


