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
CAdaBoost::CAdaBoost(SEXP radMisc): CDistribution(radMisc)
{
}

//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CAdaBoost::Create(SEXP radMisc,
		const char* szIRMeasure, int& cTrain)
{
 	return new CAdaBoost(radMisc);
}

CAdaBoost::~CAdaBoost()
{
}

void CAdaBoost::ComputeWorkingResponse
(
 const CDataset* pData,
 const double *adF,
 double *adZ
)
{


	for(long i=0; i<pData->get_trainSize(); i++)
	{
		adZ[i] = -(2*pData->y_ptr()[i]-1) * std::exp(-(2*pData->y_ptr()[i]-1)*(pData->offset_ptr(false)[i]+adF[i]));
	}


}



double CAdaBoost::InitF
(
 const CDataset* pData
)
{
    double dNum = 0.0;
    double dDen = 0.0;


	for(long i=0; i< pData->get_trainSize(); i++)
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

    
    return 0.5*std::log(dNum/dDen);
}


double CAdaBoost::Deviance
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



	for(i=0; i!=cLength; i++)
	{
		dL += pData->weight_ptr()[i] * std::exp(-(2*pData->y_ptr()[i]-1)*(pData->offset_ptr(false)[i]+adF[i]));
		dW += pData->weight_ptr()[i];
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
	const CDataset* pData,
    const double *adF,
    unsigned long cTermNodes,
    double* adZ,
    CTreeComps* pTreeComps
)
{
  double dF = 0.0;
  unsigned long iObs = 0;
  unsigned long iNode = 0;
  vecdNum.resize(cTermNodes);
  vecdNum.assign(vecdNum.size(),0.0);
  vecdDen.resize(cTermNodes);
  vecdDen.assign(vecdDen.size(),0.0);
    

  for(iObs=0; iObs< pData->get_trainSize(); iObs++)
    {
      if(pData->GetBagElem(iObs))
        {
	  dF = adF[iObs] + ((pData->offset_ptr(false)==NULL) ? 0.0 : pData->offset_ptr(false)[iObs]);
	  vecdNum[pTreeComps->GetNodeAssign()[iObs]] +=
	    pData->weight_ptr()[iObs]*(2*pData->y_ptr()[iObs]-1)*std::exp(-(2*pData->y_ptr()[iObs]-1)*dF);
	  vecdDen[pTreeComps->GetNodeAssign()[iObs]] +=
	    pData->weight_ptr()[iObs]*std::exp(-(2*pData->y_ptr()[iObs]-1)*dF);
        }
    }
  
  for(iNode=0; iNode<cTermNodes; iNode++)
    {
      if(pTreeComps->GetTermNodes()[iNode]!=NULL)
        {
	  if(vecdDen[iNode] == 0)
            {
	      	  pTreeComps->GetTermNodes()[iNode]->dPrediction = 0.0;
            }
	  else
            {
	      pTreeComps->GetTermNodes()[iNode]->dPrediction =
		vecdNum[iNode]/vecdDen[iNode];
            }
        }
    }
}


double CAdaBoost::BagImprovement
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
        if(!data.GetBag()[i])
        {
            dF = adF[i] + data.offset_ptr(false)[i];

            dReturnValue += data.weight_ptr()[i]*
                (std::exp(-(2*data.y_ptr()[i]-1)*dF) -
                 std::exp(-(2*data.y_ptr()[i]-1)*(dF+(shrinkage)*(adFadj[i]))));
            dW += data.weight_ptr()[i];
        }
    }

    return dReturnValue/dW;
}
