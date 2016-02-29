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
CGaussian::CGaussian(SEXP radMisc, const CDataset& data): CDistribution(radMisc, data)
{
}

//----------------------------------------
// Function Members - Public
//----------------------------------------
std::auto_ptr<CDistribution> CGaussian::Create(SEXP radMisc, const CDataset& data,
											const char* szIRMeasure, int& cGroups, int& cTrain)
{
	return std::auto_ptr<CDistribution>(new CGaussian(radMisc, data));
}

CGaussian::~CGaussian()
{
}


void CGaussian::ComputeWorkingResponse
(
 const double *adF,
 double *adZ,
 const bag& afInBag,
 unsigned long nTrain
 )
{
  unsigned long i = 0;
  
  if (!(pData->y_ptr() && adF && adZ && pData->weight_ptr())) {
    throw GBM::invalid_argument();
  }
  
  if(pData->offset_ptr(false) == NULL)
    {
      for(i=0; i<nTrain; i++)
        {
	  adZ[i] = pData->y_ptr()[i] - adF[i];
        }
    }
  else
    {
      for(i=0; i<nTrain; i++)
        {
	  adZ[i] = pData->y_ptr()[i] - pData->offset_ptr(false)[i] - adF[i];
        }
    }
}

void CGaussian::InitF
(
    double &dInitF,
    unsigned long cLength
)
{
    double dSum=0.0;
    double dTotalWeight = 0.0;
    unsigned long i=0;

    // compute the mean
    if(pData->offset_ptr(false)==NULL)
    {
        for(i=0; i<cLength; i++)
        {
            dSum += pData->weight_ptr()[i]*pData->y_ptr()[i];
            dTotalWeight += pData->weight_ptr()[i];
        }
    }
    else
    {
        for(i=0; i<cLength; i++)
        {
            dSum += pData->weight_ptr()[i]*(pData->y_ptr()[i] - pData->offset_ptr(false)[i]);
            dTotalWeight += pData->weight_ptr()[i];
        }
    }
    dInitF = dSum/dTotalWeight;
}


double CGaussian::Deviance
(
    const double *adF,
    unsigned long cLength,
    bool isValidationSet
)
{
    unsigned long i=0;
    double dL = 0.0;
    double dW = 0.0;

    if(isValidationSet)
    {
    	pData->shift_to_validation();
    }

    if(pData->offset_ptr(false) == NULL)
    {
        for(i=0; i<cLength; i++)
        {
            dL += pData->weight_ptr()[i]*(pData->y_ptr()[i]-adF[i])*(pData->y_ptr()[i]-adF[i]);
            dW += pData->weight_ptr()[i];
        }
    }
    else
    {
        for(i=0; i<cLength; i++)
        {
            dL += pData->weight_ptr()[i]*(pData->y_ptr()[i]-pData->offset_ptr(false)[i]-adF[i])*
                              (pData->y_ptr()[i]-pData->offset_ptr(false)[i]-adF[i]);
            dW += pData->weight_ptr()[i];
       }
    }

    if(isValidationSet)
    {
    	pData->shift_to_train();
    }

    return dL/dW;
}


void CGaussian::FitBestConstant
(
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
  // the tree aready stores the mean prediction
  // no refitting necessary
}

double CGaussian::BagImprovement
(
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

            dReturnValue += pData->weight_ptr()[i]*dStepSize*adFadj[i]*
                            (2.0*(pData->y_ptr()[i]-dF) - dStepSize*adFadj[i]);
            dW += pData->weight_ptr()[i];
        }
    }

    return dReturnValue/dW;
}



