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
CQuantile::CQuantile(SEXP radMisc, const CDataset& data): CDistribution(radMisc, data),
mpLocM("Other")
{
	dAlpha = CDistribution::misc_ptr(true)[0];
}


//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CQuantile::Create(SEXP radMisc, const CDataset& data,
											const char* szIRMeasure,
											int& cGroups, int& cTrain)
{
	return new CQuantile(radMisc, data);
}

CQuantile::~CQuantile()
{
}

void CQuantile::ComputeWorkingResponse
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
            adZ[i] = (pData->y_ptr()[i] > adF[i]) ? dAlpha : -(1.0-dAlpha);
        }
    }
    else
    {
        for(i=0; i<nTrain; i++)
        {
            adZ[i] = (pData->y_ptr()[i] > adF[i]+pData->offset_ptr(false)[i]) ? dAlpha : -(1.0-dAlpha);
        }
    }
}


void CQuantile::InitF
(
    double &dInitF,
    unsigned long cLength
)
{
    double dOffset=0.0;
    unsigned long i=0;
    int nLength = int(cLength);

    vecd.resize(cLength);
    for(i=0; i<cLength; i++)
    {
        dOffset = (pData->offset_ptr(false)==NULL) ? 0.0 : pData->offset_ptr(false)[i];
        vecd[i] = pData->y_ptr()[i] - dOffset;
    }

    dInitF = mpLocM.weightedQuantile(nLength, &vecd[0], pData->weight_ptr(), dAlpha);
}


double CQuantile::Deviance
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

    if(pData->offset_ptr(false) == NULL)
    {
        for(i=0; i<cLength; i++)
        {
            if(pData->y_ptr()[i] > adF[i])
            {
                dL += pData->weight_ptr()[i]*dAlpha*(pData->y_ptr()[i] - adF[i]);
            }
            else
            {
                dL += pData->weight_ptr()[i]*(1.0-dAlpha)*(adF[i] - pData->y_ptr()[i]);
            }
            dW += pData->weight_ptr()[i];
        }
    }
    else
    {
        for(i=0; i<cLength; i++)
        {
            if(pData->y_ptr()[i] > adF[i] + pData->offset_ptr(false)[i])
            {
                dL += pData->weight_ptr()[i]*dAlpha      *(pData->y_ptr()[i] - adF[i]-pData->offset_ptr(false)[i]);
            }
            else
            {
                dL += pData->weight_ptr()[i]*(1.0-dAlpha)*(adF[i]+pData->offset_ptr(false)[i] - pData->y_ptr()[i]);
            }
            dW += pData->weight_ptr()[i];
        }
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
    const double *adF,
    double *adZ,
    const std::vector<unsigned long> &aiNodeAssign,
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

  vecd.resize(nTrain); // should already be this size from InitF
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
		  
		  vecd[iVecd] = pData->y_ptr()[iObs] - dOffset - adF[iObs];
		  adW2[iVecd] = pData->weight_ptr()[iObs];
		  iVecd++;
                }
            }
	  
	  vecpTermNodes[iNode]->dPrediction = mpLocM.weightedQuantile(iVecd, &vecd[0], &adW2[0], dAlpha);
	}
    }
}



double CQuantile::BagImprovement
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
            if(pData->y_ptr()[i] > dF)
            {
                dReturnValue += pData->weight_ptr()[i]*dAlpha*(pData->y_ptr()[i]-dF);
            }
            else
            {
                dReturnValue += pData->weight_ptr()[i]*(1-dAlpha)*(dF-pData->y_ptr()[i]);
            }

            if(pData->y_ptr()[i] > dF+dStepSize*adFadj[i])
            {
                dReturnValue -= pData->weight_ptr()[i]*dAlpha*
                                (pData->y_ptr()[i] - dF-dStepSize*adFadj[i]);
            }
            else
            {
                dReturnValue -= pData->weight_ptr()[i]*(1-dAlpha)*
                                (dF+dStepSize*adFadj[i] - pData->y_ptr()[i]);
            }
            dW += pData->weight_ptr()[i];
        }
    }

    return dReturnValue/dW;
}

