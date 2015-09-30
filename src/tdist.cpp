//  GBM by Greg Ridgeway  Copyright (C) 2003

#include <vector>

#include "tdist.h"

void CTDist::ComputeWorkingResponse(
     const double *adY,
     const double *adMisc,
     const double *adOffset,
     const double *adF,
     double *adZ,
     const double *adWeight,
     const bag& afInBag,
     unsigned long nTrain){

  unsigned long i = 0;
  double dU = 0.0;

  if(adOffset == NULL) {
    for(i=0; i<nTrain; i++) {
      dU = adY[i] - adF[i];
      adZ[i] = (2 * dU) / (mdNu + (dU * dU));
    }
  }
  else {
    for(i=0; i<nTrain; i++) {
      dU = adY[i] - adOffset[i] - adF[i];
      adZ[i] = (2 * dU) / (mdNu + (dU * dU));
    }
  }
} // Close void CTDist::ComputeWorkingResponse

void CTDist::InitF(
     const double *adY,
     const double *adMisc,
     const double *adOffset,
     const double *adWeight,
     double &dInitF,
     unsigned long cLength) {

  // Local variables
  int ii;

  // Get objects to pass into the LocM function
  int iN = int(cLength);
  std::vector<double> adArr(cLength);

  for (ii = 0; ii < iN; ii++) {
    double dOffset = (adOffset==NULL) ? 0.0 : adOffset[ii];
    adArr[ii] = adY[ii] - dOffset;
  }

  dInitF = mpLocM.LocationM(iN, &adArr[0], adWeight, 0.5);
} // Close CTDist::InitF

double CTDist::Deviance(
       const double *adY,
       const double *adMisc,
       const double *adOffset,
       const double *adWeight,
       const double *adF,
       unsigned long cLength) {
  unsigned long i=0;
  double dL = 0.0;
  double dW = 0.0;
  double dU = 0.0;

  if(adOffset == NULL) {
    for(i=0; i<cLength; i++) {
      dU = adY[i] - adF[i];
      dL += adWeight[i] * std::log(mdNu + (dU * dU));
      dW += adWeight[i];
    }
  }
  else {
    for(i=0; i<cLength; i++) {
      dU = adY[i] - adOffset[i] - adF[i];
      dL += adWeight[i] * std::log(mdNu + (dU * dU));
      dW += adWeight[i];
    }
  }
  return dL/dW;
} // Close CTDist::Deviance

void CTDist::FitBestConstant(
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
     const double *adFadj) {
  // Local variables
  unsigned long iNode = 0;
  unsigned long iObs = 0;

  std::vector<double> adArr, adWeight;
  // Call LocM for the array of values on each node
  for(iNode=0; iNode<cTermNodes; iNode++) {
    if(vecpTermNodes[iNode] -> cN >= cMinObsInNode) {
      adArr.clear();
      adWeight.clear();

      for (iObs = 0; iObs < nTrain; iObs++) {
        if(afInBag[iObs] && (aiNodeAssign[iObs] == iNode)) {
          const double dOffset = (adOffset==NULL) ? 0.0 : adOffset[iObs];
          adArr.push_back(adY[iObs] - dOffset - adF[iObs]);
          adWeight.push_back(adW[iObs]);
        }
      }
      
      vecpTermNodes[iNode] -> dPrediction = mpLocM.LocationM(adArr.size(), &adArr[0],
                                                             &adWeight[0], 0.5);
    } // Close if (vecpTermNodes)
  } // Close for(iNode=0)
} // Close CTDist:FitBestConstant

double CTDist::BagImprovement(
       const double *adY,
       const double *adMisc,
       const double *adOffset,
       const double *adWeight,
       const double *adF,
       const double *adFadj,
       const bag& afInBag,
       double dStepSize,
       unsigned long nTrain){
  double dReturnValue = 0.0;
  unsigned long i = 0;
  double dW = 0.0;

  for(i=0; i<nTrain; i++) {
    if(!afInBag[i]) {
      const double dF = adF[i] + ((adOffset==NULL) ? 0.0 : adOffset[i]);
      const double dU = (adY[i] - dF);
      const double dV = (adY[i] - dF - dStepSize * adFadj[i]) ;

      dReturnValue += adWeight[i] * (std::log(mdNu + (dU * dU)) - log(mdNu + (dV * dV)));
      dW += adWeight[i];
    }
  }

    return dReturnValue/dW;
}
