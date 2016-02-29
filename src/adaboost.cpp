//-----------------------------------
//
// File: adaboost.cpp
//
// Description: distribution used for adaboosting.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "adaboost.h"
#include <memory>

//----------------------------------------
// Function Members - Private
//----------------------------------------
CAdaBoost::CAdaBoost(SEXP radMisc, const CDataset& data): CDistribution(radMisc, data)
{
}

//----------------------------------------
// Function Members - Public
//----------------------------------------
std::auto_ptr<CDistribution> CAdaBoost::Create(SEXP radMisc, const CDataset& data,
		const char* szIRMeasure,
		int& cGroups, int& cTrain)
{
 	return std::auto_ptr<CDistribution>(new CAdaBoost(radMisc, data));
}

CAdaBoost::~CAdaBoost()
{
}

void CAdaBoost::ComputeWorkingResponse
(
 const double *adF,
 double *adZ,
 const bag& afInBag,
 unsigned long nTrain
)
{
  unsigned long i = 0;
  if(pData->offset_ptr(false) == NULL)
    {
      for(i=0; i<nTrain; i++)
        {
	  adZ[i] = -(2*pData->y_ptr()[i]-1) * std::exp(-(2*pData->y_ptr()[i]-1)*adF[i]);
        }
    }
  else
    {
      for(i=0; i<nTrain; i++)
        {
	  adZ[i] = -(2*pData->y_ptr()[i]-1) * std::exp(-(2*pData->y_ptr()[i]-1)*(pData->offset_ptr(false)[i]+adF[i]));
        }
    }

}



void CAdaBoost::InitF
(
 double &dInitF,
 unsigned long cLength
)
{
    unsigned long i=0;
    double dNum = 0.0;
    double dDen = 0.0;

    dInitF = 0.0;

    if(pData->offset_ptr(false) == NULL)
    {
      for(i=0; i<cLength; i++)
        {
	  if(pData->y_ptr()[i]==1.0)
            {
	      dNum += pData->weight_ptr()[i];
            }
	  else
            {
                dDen += pData->weight_ptr()[i];
            }
        }
    }
    else
      {
        for(i=0; i<cLength; i++)
        {
	  if(pData->y_ptr()[i]==1.0)
            {
	      dNum += pData->weight_ptr()[i] * std::exp(-pData->offset_ptr(false)[i]);
            }
	  else
            {
	      dDen += pData->weight_ptr()[i] * std::exp(pData->offset_ptr(false)[i]);
            }
        }
      }
    
    dInitF = 0.5*std::log(dNum/dDen);
}


double CAdaBoost::Deviance
(
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
    else
    {
    	const double* adY = pData->y_ptr();
    	const double* adW = pData->weight_ptr();
    	const double* adOff = pData->offset_ptr(false);
    }

    if(pData->offset_ptr(false) == NULL)
    {
        for(i=0; i!=cLength; i++)
        {
            dL += pData->weight_ptr()[i] * std::exp(-(2*pData->y_ptr()[i]-1)*adF[i]);
            dW += pData->weight_ptr()[i];
        }
    }
    else
    {
        for(i=0; i!=cLength; i++)
        {
            dL += pData->weight_ptr()[i] * std::exp(-(2*pData->y_ptr()[i]-1)*(pData->offset_ptr(false)[i]+adF[i]));
            dW += pData->weight_ptr()[i];
       }
    }

    // Switch back to trainig set if necessary
   if(isValidationSet)
   {
	   pData->shift_to_train();
   }

    return dL/dW;
}


void CAdaBoost::FitBestConstant
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
	  dF = adF[iObs] + ((pData->offset_ptr(false)==NULL) ? 0.0 : pData->offset_ptr(false)[iObs]);
	  vecdNum[aiNodeAssign[iObs]] +=
	    pData->weight_ptr()[iObs]*(2*pData->y_ptr()[iObs]-1)*std::exp(-(2*pData->y_ptr()[iObs]-1)*dF);
	  vecdDen[aiNodeAssign[iObs]] +=
	    pData->weight_ptr()[iObs]*std::exp(-(2*pData->y_ptr()[iObs]-1)*dF);
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
                (std::exp(-(2*pData->y_ptr()[i]-1)*dF) -
                 std::exp(-(2*pData->y_ptr()[i]-1)*(dF+dStepSize*adFadj[i])));
            dW += pData->weight_ptr()[i];
        }
    }

    return dReturnValue/dW;
}
