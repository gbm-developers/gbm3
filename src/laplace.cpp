//  GBM by Greg Ridgeway  Copyright (C) 2003

#include <vector>
#include "laplace.h"


CLaplace::~CLaplace()
{
}


void CLaplace::ComputeWorkingResponse
(
 const double *adY,
 const double *adMisc,
 const double *adOffset,
 const double *adF,
 double *adZ,
 const double *adWeight,
 const bag& afInBag,
 unsigned long nTrain,
 int cIdxOff
)
{
    unsigned long i = 0;

    if(adOffset == NULL)
    {
        for(i=0; i<nTrain; i++)
        {
            adZ[i] = (adY[i] - adF[i]) > 0.0 ? 1.0 : -1.0;
        }
    }
    else
    {
        for(i=0; i<nTrain; i++)
        {
            adZ[i] = (adY[i] - adOffset[i] - adF[i]) > 0.0 ? 1.0 : -1.0;
        }
    }
}



void CLaplace::InitF
(
    const double *adY,
    const double *adMisc,
    const double *adOffset,
    const double *adWeight,
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
      dOffset = (adOffset==NULL) ? 0.0 : adOffset[ii];
      adArr[ii] = adY[ii] - dOffset;
    }
  
  dInitF = mpLocM.weightedQuantile(nLength, &adArr[0], adWeight, 0.5); // median
}


double CLaplace::Deviance
(
    const double *adY,
    const double *adMisc,
    const double *adOffset,
    const double *adWeight,
    const double *adF,
    unsigned long cLength,
   int cIdxOff
)
{
    unsigned long i=0;
    double dL = 0.0;
    double dW = 0.0;

    if(adOffset == NULL)
    {
        for(i=cIdxOff; i<cLength+cIdxOff; i++)
        {
            dL += adWeight[i]*fabs(adY[i]-adF[i]);
            dW += adWeight[i];
        }
    }
    else
    {
        for(i=cIdxOff; i<cLength+cIdxOff; i++)
        {
            dL += adWeight[i]*fabs(adY[i]-adOffset[i]-adF[i]);
            dW += adWeight[i];
        }
    }

    return dL/dW;
}


// DEBUG: needs weighted median
void CLaplace::FitBestConstant
(
 const double *adY,
 const double *adMisc,
 const double *adOffset,
 const double *adW,
 const double *adF,
 double *adZ,
 const std::vector<unsigned long>& aiNodeAssign,
 unsigned long nTrain,
 VEC_P_NODETERMINAL vecpTermNodes,
 unsigned long cTermNodes,
 unsigned long cMinObsInNode,
 const bag& afInBag,
 const double *adFadj,
 int cIdxOff
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
		  dOffset = (adOffset==NULL) ? 0.0 : adOffset[iObs];
		  adArr[iVecd] = adY[iObs] - dOffset - adF[iObs];
		  adW2[iVecd] = adW[iObs];
		  iVecd++;
		}
	      
            }
	  
	  vecpTermNodes[iNode]->dPrediction = mpLocM.weightedQuantile(iVecd, &adArr[0], &adW2[0], 0.5); // median
	  
        }
    }
}



double CLaplace::BagImprovement
(
    const double *adY,
    const double *adMisc,
    const double *adOffset,
    const double *adWeight,
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
            dF = adF[i] + ((adOffset==NULL) ? 0.0 : adOffset[i]);

            dReturnValue +=
                adWeight[i]*(fabs(adY[i]-dF) - fabs(adY[i]-dF-dStepSize*adFadj[i]));
            dW += adWeight[i];
        }
    }

    return dReturnValue/dW;
}
