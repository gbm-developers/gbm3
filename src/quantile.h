//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//  File:       quantile.h
//
//  License:    GNU GPL (version 2 or later)
//
//  Contents:   laplace object
//
//  Owner:      gregr@rand.org
//
//  History:    10/8/2006   Created by Brian Kriegler (bk@stat.ucla.edu)
//              6/11/2007   gregr merged with official gbm
//
//------------------------------------------------------------------------------

#ifndef QUANTILE_H
#define QUANTILE_H

#include <algorithm>
#include "distribution.h"
#include "locationm.h"


class CQuantile: public CDistribution
{

public:

 CQuantile(double dAlpha) : dAlpha(dAlpha), mpLocM("Other") {};

  virtual ~CQuantile() {};
  
  
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
		       const std::vector<unsigned long>& aiNodeAssign,
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
    vector<double> vecd;
    double dAlpha;
    CLocationM mpLocM;
};

#endif // QUANTILE_H



