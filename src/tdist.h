//------------------------------------------------------------------------------
//  GBM alteration by Daniel Edwards
//
//  File:       tdist.h
//
//  Contains:   Distribution object to implement t-distribution
//
//  History:    04/04/2008   Created
//
//------------------------------------------------------------------------------

#ifndef TDISTCGBM_H
#define TDISTCGBM_H

#include <algorithm>
#include "distribution.h"
#include "locationm.h"


class CTDist : public CDistribution
{
  
 public:
  
 CTDist(double adNu) : mdNu(adNu), mpLocM("tdist", adNu) {};

  virtual ~CTDist() {};
  
  
  void ComputeWorkingResponse(double *adY,
			      double *adMisc,
			      double *adOffset,
			      double *adF,
			      double *adZ,
			      double *adWeight,
			      int *afInBag,
			      unsigned long nTrain,
			      int cIdxOff);

    void InitF(double *adY,
	       double *adMisc,
	       double *adOffset,
	       double *adWeight,
	       double &dInitF,
	       unsigned long cLength);
    
    void FitBestConstant(double *adY,
			 double *adMisc,
			 double *adOffset,
			 double *adW,
			 double *adF,
			 double *adZ,
			 const std::vector<unsigned long> &aiNodeAssign,
			 unsigned long nTrain,
			 VEC_P_NODETERMINAL vecpTermNodes,
			 unsigned long cTermNodes,
			 unsigned long cMinObsInNode,
			 int *afInBag,
			 double *adFadj,
			 int cIdxOff);

    double Deviance(double *adY,
                    double *adMisc,
                    double *adOffset,
                    double *adWeight,
                    double *adF,
                    unsigned long cLength,
		    int cIdxOff);

    double BagImprovement(double *adY,
                          double *adMisc,
                          double *adOffset,
                          double *adWeight,
                          double *adF,
                          double *adFadj,
                          int *afInBag,
                          double dStepSize,
                          unsigned long nTrain);

private:
    double mdNu;
    CLocationM mpLocM;
};

#endif // TDISTCGBM_H



