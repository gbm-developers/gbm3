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
  
  
  void ComputeWorkingResponse(const double *adY,
			      const double *adMisc,
			      const double *adOffset,
			      const double *adF,
			      double *adZ,
			      const double *adWeight,
			      const bag& afInBag,
			      unsigned long nTrain);

    void InitF(const double *adY,
	       const double *adMisc,
	       const double *adOffset,
	       const double *adWeight,
	       double &dInitF,
	       unsigned long cLength);
    
    void FitBestConstant(const double *adY,
			 const double *adMisc,
			 const double *adOffset,
			 const double *adW,
			 const double *adF,
			 double *adZ,
			 const std::vector<unsigned long> &aiNodeAssign,
			 unsigned long nTrain,
			 VEC_P_NODETERMINAL vecpTermNodes,
			 unsigned long cTermNodes,
			 unsigned long cMinObsInNode,
			 const bag& afInBag,
			 const double *adFadj);

    double Deviance(const double *adY,
                    const double *adMisc,
                    const double *adOffset,
                    const double *adWeight,
                    const double *adF,
                    unsigned long cLength);

    double BagImprovement(const double *adY,
                          const double *adMisc,
                          const double *adOffset,
                          const double *adWeight,
                          const double *adF,
                          const double *adFadj,
                          const bag& afInBag,
                          double dStepSize,
                          unsigned long nTrain);

private:
    double mdNu;
    CLocationM mpLocM;
};

#endif // TDISTCGBM_H



