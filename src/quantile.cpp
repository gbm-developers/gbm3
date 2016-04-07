//-----------------------------------
//
// File: quantile.cpp
//
// Description: quantile distribution for GBM.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "quantile.h"

//----------------------------------------
// Function Members - Private
//----------------------------------------
CQuantile::CQuantile(SEXP radMisc): CDistribution(radMisc),
mpLocM("Other")
{
	dAlpha = CDistribution::misc_ptr(true)[0];
}


//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CQuantile::Create(SEXP radMisc,
								 const char* szIRMeasure,
								int& cTrain)
{
	return new CQuantile(radMisc);
}

CQuantile::~CQuantile()
{
}

void CQuantile::ComputeWorkingResponse
(
	const CDataset* pData,
    const double *adF,
    double *adZ
)
{
    unsigned long i = 0;
	for(i=0; i<pData->get_trainSize(); i++)
	{
		adZ[i] = (pData->y_ptr()[i] > adF[i]+pData->offset_ptr(false)[i]) ? dAlpha : -(1.0-dAlpha);
	}

}


double CQuantile::InitF
(
	const CDataset* pData
)
{
    double dOffset=0.0;
    vecd.resize(pData->get_trainSize());
    for(long i=0; i< pData->get_trainSize(); i++)
    {
        dOffset = pData->offset_ptr(false)[i];
        vecd[i] = pData->y_ptr()[i] - dOffset;
    }

    return mpLocM.weightedQuantile(pData->get_trainSize(), &vecd[0], pData->weight_ptr(), dAlpha);
}


double CQuantile::Deviance
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
		if(pData->y_ptr()[i] > adF[i] + pData->offset_ptr(false)[i])
		{
			dL += pData->weight_ptr()[i]*dAlpha*(pData->y_ptr()[i] - adF[i]-pData->offset_ptr(false)[i]);
		}
		else
		{
			dL += pData->weight_ptr()[i]*(1.0-dAlpha)*(adF[i]+pData->offset_ptr(false)[i] - pData->y_ptr()[i]);
		}
		dW += pData->weight_ptr()[i];
	}


    // Switch back to training set if necessary
    if(isValidationSet)
    {
 	   pData->shift_to_train();
    }

    return dL/dW;
}

void CQuantile::FitBestConstant
(
	const CDataset* pData,
    const double *adF,
    unsigned long cTermNodes,
    double* adZ,
    CTreeComps* pTreeComps
)
{
  unsigned long iNode = 0;
  unsigned long iObs = 0;
  unsigned long iVecd = 0;
  double dOffset;

  vecd.resize(pData->get_trainSize()); // should already be this size from InitF
  std::vector<double> adW2(pData->get_trainSize());

  for(iNode=0; iNode<cTermNodes; iNode++)
    {
      if(pTreeComps->GetTermNodes()[iNode]->cN >= pTreeComps->GetMinNodeObs())
        {
	  iVecd = 0;
	  for(iObs=0; iObs< pData->get_trainSize(); iObs++)
            {
	      if(pData->GetBagElem(iObs) && (pTreeComps->GetNodeAssign()[iObs] == iNode))
                {
		  dOffset = (pData->offset_ptr(false)==NULL) ? 0.0 : pData->offset_ptr(false)[iObs];
		  
		  vecd[iVecd] = pData->y_ptr()[iObs] - dOffset - adF[iObs];
		  adW2[iVecd] = pData->weight_ptr()[iObs];
		  iVecd++;
                }
            }
	  
	 pTreeComps->GetTermNodes()[iNode]->dPrediction = mpLocM.weightedQuantile(iVecd, &vecd[0], &adW2[0], dAlpha);
	}
    }
}


double CQuantile::BagImprovement
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

    for(i=0; i <data.get_trainSize(); i++)
    {
        if(!data.GetBagElem(i))
        {
            dF = adF[i] + data.offset_ptr(false)[i];

            if(data.y_ptr()[i] > dF)
            {
                dReturnValue += data.weight_ptr()[i]*dAlpha*(data.y_ptr()[i]-dF);
            }
            else
            {
                dReturnValue += data.weight_ptr()[i]*(1-dAlpha)*(dF-data.y_ptr()[i]);
            }

            if(data.y_ptr()[i] > dF+shrinkage*adFadj[i])
            {
                dReturnValue -= data.weight_ptr()[i]*dAlpha*
                                (data.y_ptr()[i] - dF-shrinkage*adFadj[i]);
            }
            else
            {
                dReturnValue -= data.weight_ptr()[i]*(1-dAlpha)*
                                (dF+shrinkage*adFadj[i] - data.y_ptr()[i]);
            }
            dW += data.weight_ptr()[i];
        }
    }
    return dReturnValue/dW;
}

