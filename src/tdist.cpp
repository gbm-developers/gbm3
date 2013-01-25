//  GBM by Greg Ridgeway  Copyright (C) 2003

#include "tdist.h"

CTDist::CTDist(double adNu)
{
    mdNu = adNu;

	double *adParams = new double[1];
	adParams[0] = adNu;

	mpLocM = new CLocationM("tdist", 1, adParams);

	delete[] adParams;
}

CTDist::~CTDist()
{
	delete mpLocM;
}


GBMRESULT CTDist::ComputeWorkingResponse
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adF, 
    double *adZ, 
    double *adWeight,
    bool *afInBag,
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

    return GBM_OK;
}


GBMRESULT CTDist::InitF
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
	double *adArr = new double[iN];

	for (ii = 0; ii < iN; ii++)
	{
		double dOffset = (adOffset==NULL) ? 0.0 : adOffset[ii];
		adArr[ii] = adY[ii] - dOffset;
	}

	dInitF = mpLocM->LocationM(iN, adArr, adWeight);

    delete[] adArr;

    return GBM_OK;
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
			dL += adWeight[i] * log(mdNu + (dU * dU));
            dW += adWeight[i];
        }
    }
    else
    {
        for(i=cIdxOff; i<cLength+cIdxOff; i++)
        {
			dU = adY[i] - adOffset[i] - adF[i];
		    dL += adWeight[i] * log(mdNu + (dU * dU));
            dW += adWeight[i];
        }
    }

    return dL/dW;
}


GBMRESULT CTDist::FitBestConstant
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adW,
    double *adF,
    double *adZ,
    unsigned long *aiNodeAssign,
    unsigned long nTrain,
    VEC_P_NODETERMINAL vecpTermNodes,
    unsigned long cTermNodes,
    unsigned long cMinObsInNode,
    bool *afInBag,
    double *adFadj,
	int cIdxOff
)
{
   	// Local variables
    GBMRESULT hr = GBM_OK;
	unsigned long iNode = 0;
    unsigned long iObs = 0;
    double dOffset;

	// Call LocM for the array of values on each node
    for(iNode=0; iNode<cTermNodes; iNode++)
    {
        if(vecpTermNodes[iNode]->cN >= cMinObsInNode)
        {
			// Get the number of nodes here
			int iNumNodes = 0;
			for (iObs = 0; iObs < nTrain; iObs++)
			{
				if(afInBag[iObs] && (aiNodeAssign[iObs] == iNode))
                {
                    iNumNodes++;
                }
			}

			// Create the arrays to centre
			double *adArr = new double[iNumNodes];
			double *adWeight = new double[iNumNodes];

			int iIdx = 0;
			for(iObs=0; iObs<nTrain; iObs++)
            {
                if(afInBag[iObs] && (aiNodeAssign[iObs] == iNode))
                {
                    dOffset = (adOffset==NULL) ? 0.0 : adOffset[iObs];
                    adArr[iIdx] = adY[iObs] - dOffset - adF[iObs];
					adWeight[iIdx] = adW[iObs];

					iIdx++;
                }
            }

           	vecpTermNodes[iNode]->dPrediction = mpLocM->LocationM(iNumNodes, adArr, 
				                                                 adWeight); 

			delete[] adArr;
			delete[] adWeight;

        }
    }

    return hr;
}

double CTDist::BagImprovement
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adWeight,
    double *adF,
    double *adFadj,
    bool *afInBag,
    double dStepSize,
    unsigned long nTrain
)
{
    double dReturnValue = 0.0;
    double dF = 0.0;
    double dW = 0.0;
    unsigned long i = 0;
	double dU = 0.0;
    double dV = 0.0;

    for(i=0; i<nTrain; i++)
    {
        if(!afInBag[i])
        {
            dF = adF[i] + ((adOffset==NULL) ? 0.0 : adOffset[i]);
            
			dU = (adY[i] - dF);
			dV = (adY[i] - dF - dStepSize * adFadj[i]) ;

            dReturnValue += adWeight[i] * (log(mdNu + (dU * dU)) - log(mdNu + (dV * dV)));
            dW += adWeight[i];
        }
    }

    return dReturnValue/dW;
}




