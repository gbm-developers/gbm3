//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//  File:       laplace.h
//
//  License:    GNU GPL (version 2 or later)
//
//  Contents:   laplace object
//
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef LAPLACGBM_H
#define LAPLACGBM_H

#include <algorithm>
#include "distribution.h"
#include "locationm.h"


class CLaplace : public CDistribution
{

public:

 CLaplace()  : mpLocM("Other") {};

  virtual ~CLaplace();
  
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
  vector<double> vecd;
  vector<double>::iterator itMedian;
  CLocationM mpLocM;
};

#endif // LAPLACGBM_H



