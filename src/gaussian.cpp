//-----------------------------------
//
// File: gaussian.cpp
//
// Description: gaussian distribution implementation for GBM.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "gaussian.h"

//----------------------------------------
// Function Members - Private
//----------------------------------------
CGaussian::CGaussian(SEXP radMisc): CDistribution(radMisc)
{

}

//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CGaussian::Create(SEXP radMisc,
											const char* szIRMeasure, int& cTrain)
{
	return new CGaussian(radMisc);
}

CGaussian::~CGaussian()
{
}


void CGaussian::ComputeWorkingResponse
(
 const CDataset* pData,
 const double *adF,
 double *adZ
 )
{
  unsigned long i = 0;
  
  if (!(pData->y_ptr() && adF && adZ && pData->weight_ptr())) {
    throw GBM::invalid_argument();
  }
  

	for(i=0; i<pData->get_trainSize(); i++)
	{
		adZ[i] = pData->y_ptr()[i] - pData->offset_ptr(false)[i] - adF[i];
	}

}

double CGaussian::InitF
(
	const CDataset* pData
)
{
    double dSum=0.0;
    double dTotalWeight = 0.0;
    unsigned long i=0;

    // compute the mean

	for(i=0; i<pData->get_trainSize(); i++)
	{
		dSum += pData->weight_ptr()[i]*(pData->y_ptr()[i] - pData->offset_ptr(false)[i]);
		dTotalWeight += pData->weight_ptr()[i];
	}


    return dSum/dTotalWeight;
}


double CGaussian::Deviance
(
	const CDataset* pData,
    const double *adF,
    bool isValidationSet
)
{
    unsigned long i=0;
    double dL = 0.0;
    double dW = 0.0;

    long cLength = pData->get_trainSize();
    if(isValidationSet)
    {
    	pData->shift_to_validation();
    	cLength = pData->GetValidSize();
    }



	for(i=0; i<cLength; i++)
	{
		dL += pData->weight_ptr()[i]*(pData->y_ptr()[i]-pData->offset_ptr(false)[i]-adF[i])*
						  (pData->y_ptr()[i]-pData->offset_ptr(false)[i]-adF[i]);
		dW += pData->weight_ptr()[i];
	}


    if(isValidationSet)
    {
    	pData->shift_to_train();
    }

    return dL/dW;
}


void CGaussian::FitBestConstant
(
	const CDataset* pData,
    const double *adF,
    unsigned long cTermNodes,
    double* adZ,
    CTreeComps* pTreeComps
)
{
  // the tree aready stores the mean prediction
  // no refitting necessary
}

double CGaussian::BagImprovement
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

            dReturnValue += data.weight_ptr()[i]*shrinkage*adFadj[i]*
                            (2.0*(data.y_ptr()[i]-dF) - shrinkage*adFadj[i]);
            dW += data.weight_ptr()[i];
        }
    }

    return dReturnValue/dW;
}



