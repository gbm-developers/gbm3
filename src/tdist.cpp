//  GBM by Greg Ridgeway  Copyright (C) 2003

#include <vector>

#include "tdist.h"

void CTDist::ComputeWorkingResponse
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adF,
    double *adZ,
    double *adWeight,
    int *afInBag,
    unsigned long nTrain,
    int cIdxOff
)
{
    unsigned long i = 0;
    double dU = 0.0;

    if(adOffset == NULL)
    {
        for(i=0; i<nTrain; i++)
        {
	  dU = adY[i] - adF[i];
	  adZ[i] = (2 * dU) / (mdNu + (dU * dU));
        }
    }
    else
    {
        for(i=0; i<nTrain; i++)
        {
 		    dU = adY[i] - adOffset[i] - adF[i];
			adZ[i] = (2 * dU) / (mdNu + (dU * dU));
        }
    }
}


void CTDist::InitF
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adWeight,
    double &dInitF,
    unsigned long cLength
)
{

	// Local variables
	int ii;

	// Get objects to pass into the LocM function
	int iN = int(cLength);
	std::vector<double> adArr(cLength);

	for (ii = 0; ii < iN; ii++)
	{
		double dOffset = (adOffset==NULL) ? 0.0 : adOffset[ii];
		adArr[ii] = adY[ii] - dOffset;
	}

	dInitF = mpLocM.LocationM(iN, &adArr[0], adWeight, 0.5);
}

double CTDist::Deviance
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adWeight,
    double *adF,
    unsigned long cLength,
    int cIdxOff
)
{
    unsigned long i=0;
    double dL = 0.0;
    double dW = 0.0;
	double dU = 0.0;

    if(adOffset == NULL)
    {
        for(i=cIdxOff; i<cLength+cIdxOff; i++)
        {
			dU = adY[i] - adF[i];
			dL += adWeight[i] * std::log(mdNu + (dU * dU));
            dW += adWeight[i];
        }
    }
    else
    {
        for(i=cIdxOff; i<cLength+cIdxOff; i++)
        {
			dU = adY[i] - adOffset[i] - adF[i];
		    dL += adWeight[i] * std::log(mdNu + (dU * dU));
            dW += adWeight[i];
        }
    }

    return dL/dW;
}


void CTDist::FitBestConstant
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adW,
    double *adF,
    double *adZ,
    const std::vector<unsigned long>& aiNodeAssign,
    unsigned long nTrain,
    VEC_P_NODETERMINAL vecpTermNodes,
    unsigned long cTermNodes,
    unsigned long cMinObsInNode,
    int *afInBag,
    double *adFadj,
	int cIdxOff
)
{
   	// Local variables
    unsigned long iNode = 0;
    unsigned long iObs = 0;


    std::vector<double> adArr, adWeight;
	// Call LocM for the array of values on each node
    for(iNode=0; iNode<cTermNodes; iNode++)
    {
      if(vecpTermNodes[iNode]->cN >= cMinObsInNode)
        {
	  adArr.clear();
	  adWeight.clear();

	  for (iObs = 0; iObs < nTrain; iObs++)
	    {
	      if(afInBag[iObs] && (aiNodeAssign[iObs] == iNode))
                {
		  const double dOffset = (adOffset==NULL) ? 0.0 : adOffset[iObs];
		  adArr.push_back(adY[iObs] - dOffset - adF[iObs]);
		  adWeight.push_back(adW[iObs]);
                }
	    }

	  vecpTermNodes[iNode]->dPrediction = mpLocM.LocationM(adArr.size(), &adArr[0],
							       &adWeight[0], 0.5);

        }
    }
}

double CTDist::BagImprovement
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adWeight,
    double *adF,
    double *adFadj,
    int *afInBag,
    double dStepSize,
    unsigned long nTrain
)
{
    double dReturnValue = 0.0;
    unsigned long i = 0;
    double dU = 0.0;
    double dW = 0.0;

    for(i=0; i<nTrain; i++)
    {
        if(!afInBag[i])
        {
            const double dF = adF[i] + ((adOffset==NULL) ? 0.0 : adOffset[i]);
	    const double dU = (adY[i] - dF);
	    const double dV = (adY[i] - dF - dStepSize * adFadj[i]) ;

            dReturnValue += adWeight[i] * (std::log(mdNu + (dU * dU)) - log(mdNu + (dV * dV)));
            dW += adWeight[i];
        }
    }

    return dReturnValue/dW;
}




