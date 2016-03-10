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
										int& cGroups, int& cTrain)
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
            adZ[i] = (pData->y_ptr()[i] - adF[i]) > 0.0 ? 1.0 : -1.0;
        }
    }
    else
    {
        for(i=0; i<nTrain; i++)
        {
            adZ[i] = (pData->y_ptr()[i] - pData->offset_ptr(false)[i] - adF[i]) > 0.0 ? 1.0 : -1.0;
        }
    }
}



void CLaplace::InitF
(
	const CDataset* pData,
    double &dInitF,
    unsigned long cLength
)
{
  double dOffset = 0.0;
  unsigned long ii = 0;
  int nLength = int(cLength);
  
  std::vector<double> adArr(cLength);
  
  for (ii = 0; ii < cLength; ii++)
    {
      dOffset = (pData->offset_ptr(false)==NULL) ? 0.0 : pData->offset_ptr(false)[ii];
      adArr[ii] = pData->y_ptr()[ii] - dOffset;
    }
  
  dInitF = mpLocM.weightedQuantile(nLength, &adArr[0], pData->weight_ptr(), 0.5); // median
}


double CLaplace::Deviance
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

    if(isValidationSet)
    {
    	pData->shift_to_validation();
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
  unsigned long iNode = 0;
  unsigned long iObs = 0;
  unsigned long iVecd = 0;
  double dOffset;
  
//    vecd.resize(nTrain); // should already be this size from InitF

  std::vector<double> adArr(nTrain);
  std::vector<double> adW2(nTrain);
  
  for(iNode=0; iNode<cTermNodes; iNode++)
    {
      if(vecpTermNodes[iNode]->cN >= cMinObsInNode)
        {
	  iVecd = 0;
	  for(iObs=0; iObs<nTrain; iObs++)
            {
	      if(afInBag[iObs] && (aiNodeAssign[iObs] == iNode))
                {
		  dOffset = (pData->offset_ptr(false)==NULL) ? 0.0 : pData->offset_ptr(false)[iObs];
		  adArr[iVecd] = pData->y_ptr()[iObs] - dOffset - adF[iObs];
		  adW2[iVecd] = pData->weight_ptr()[iObs];
		  iVecd++;
		}
	      
            }
	  
	  vecpTermNodes[iNode]->dPrediction = mpLocM.weightedQuantile(iVecd, &adArr[0], &adW2[0], 0.5); // median
	  
        }
    }
}



double CLaplace::BagImprovement
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

            dReturnValue +=
                pData->weight_ptr()[i]*(fabs(pData->y_ptr()[i]-dF) - fabs(pData->y_ptr()[i]-dF-dStepSize*adFadj[i]));
            dW += pData->weight_ptr()[i];
        }
    }

    return dReturnValue/dW;
}
