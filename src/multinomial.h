//------------------------------------------------------------------------------
//  GBM alteration by Daniel Edwards
//
//  File:       multinomial.h
//
//
//  Contains:   Distribution object to implement multinomial
//
//  History:    04/04/2008   Created
//
//------------------------------------------------------------------------------

#ifndef KMULTICGBM_H
#define KMULTICGBM_H

#include <algorithm>
#include <vector>
#include "distribution.h"
#include "locationm.h"

class CMultinomial : public CDistribution
{

public:

 CMultinomial(int cNumClasses, int cRows) : mcNumClasses(cNumClasses),
    mcRows(cRows), madProb(cNumClasses * cRows, 0) {};
  
  void UpdateParams(const double *adF,
		    const double *adOffset,
		    const double *adWeight,
		    unsigned long cLength);
  
  void ComputeWorkingResponse(const double *adY,
			      const double *adMisc,
			      const double *adOffset,
			      const double *adF,
			      double *adZ,
			      const double *adWeight,
			      const bag& afInBag,
			      unsigned long nTrain,
			      int cIdxOff);

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
			 const std::vector<unsigned long>& aiNodeAssign,
			 unsigned long nTrain,
			 VEC_P_NODETERMINAL vecpTermNodes,
			 unsigned long cTermNodes,
			 unsigned long cMinObsInNode,
			 const bag& afInBag,
			 const double *adFadj,
			 int cIdxOff);
    
    double Deviance(const double *adY,
                    const double *adMisc,
                    const double *adOffset,
                    const double *adWeight,
                    const double *adF,
                    unsigned long cLength,
                    int cIdxOff);

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
   unsigned long mcNumClasses;
   unsigned long mcRows;
   std::vector<double> madProb;
};

#endif // KMULTICGBM_H


