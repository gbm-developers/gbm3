//-----------------------------------
//
// File: laplace.cpp
//
// Description: laplace distribution for GBM.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "laplace.h"
#include <vector>

//----------------------------------------
// Function Members - Private
//----------------------------------------
CLaplace::CLaplace(SEXP radMisc): CDistribution(radMisc), mpLocM("Other")
{
}

//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CLaplace::Create(SEXP radMisc,
								const char* szIRMeasure,
								int& cTrain)
{
	return new CLaplace(radMisc);
}

CLaplace::~CLaplace()
{
}


void CLaplace::ComputeWorkingResponse
(
 const CDataset* pData,
 const double *adF,
 double *adZ
)
{
    unsigned long i = 0;
	for(i=0; i<pData->get_trainSize(); i++)
	{
		adZ[i] = (pData->y_ptr()[i] - pData->offset_ptr(false)[i] - adF[i]) > 0.0 ? 1.0 : -1.0;
	}

}



double CLaplace::InitF
(
	const CDataset* pData
)
{
  double dOffset = 0.0;
  unsigned long ii = 0;

  std::vector<double> adArr(pData->get_trainSize());
  
  for (ii = 0; ii < pData->get_trainSize(); ii++)
    {
      dOffset = pData->offset_ptr(false)[ii];
      adArr[ii] = pData->y_ptr()[ii] - dOffset;
    }
  
  return mpLocM.weightedQuantile(pData->get_trainSize(), &adArr[0], pData->weight_ptr(), 0.5); // median
}


double CLaplace::Deviance
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

    if(pData->offset_ptr(false) == NULL)
    {
        for(i=0; i<cLength; i++)
        {
            dL += pData->weight_ptr()[i]*fabs(pData->y_ptr()[i]-adF[i]);
            dW += pData->weight_ptr()[i];
        }
    }
    else
    {
        for(i=0; i<cLength; i++)
        {
            dL += pData->weight_ptr()[i]*fabs(pData->y_ptr()[i]-pData->offset_ptr(false)[i]-adF[i]);
            dW += pData->weight_ptr()[i];
        }
    }

    if(isValidationSet)
    {
    	pData->shift_to_train();
    }

    return dL/dW;
}


// DEBUG: needs weighted median
void CLaplace::FitBestConstant
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
  
//    vecd.resize(nTrain); // should already be this size from InitF

  std::vector<double> adArr(pData->get_trainSize());
  std::vector<double> adW2(pData->get_trainSize());
  
  for(iNode=0; iNode<cTermNodes; iNode++)
    {
      if(pTreeComps->GetTermNodes()[iNode]->cN >= pTreeComps->GetMinNodeObs())
        {
	  iVecd = 0;
	  for(iObs=0; iObs<pData->get_trainSize(); iObs++)
            {
	      if(pData->GetBagElem(iObs) && (pTreeComps->GetNodeAssign()[iObs] == iNode))
                {
		  dOffset =  pData->offset_ptr(false)[iObs];
		  adArr[iVecd] = pData->y_ptr()[iObs] - dOffset - adF[iObs];
		  adW2[iVecd] = pData->weight_ptr()[iObs];
		  iVecd++;
		}
	      
            }
	  
	  pTreeComps->GetTermNodes()[iNode]->dPrediction = mpLocM.weightedQuantile(iVecd, &adArr[0], &adW2[0], 0.5); // median
	  
        }
    }
}



double CLaplace::BagImprovement
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

            dReturnValue +=
                data.weight_ptr()[i]*(fabs(data.y_ptr()[i]-dF) - fabs(data.y_ptr()[i]-dF-shrinkage*adFadj[i]));
            dW += data.weight_ptr()[i];
        }
    }

    return dReturnValue/dW;
}
