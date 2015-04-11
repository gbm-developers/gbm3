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
  
  
  void ComputeWorkingResponse(double *adY,
			      double *adMisc,
			      double *adOffset,
			      double *adF,
			      double *adZ,
			      double *adWeight,
			      const bag& afInBag,
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
		       const std::vector<unsigned long>& aiNodeAssign,
		       unsigned long nTrain,
		       VEC_P_NODETERMINAL vecpTermNodes,
		       unsigned long cTermNodes,
		       unsigned long cMinObsInNode,
		       const bag& afInBag,
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
			const bag& afInBag,
			double dStepSize,
			unsigned long nTrain);

private:
    vector<double> vecd;
    double dAlpha;
    CLocationM mpLocM;
};

#endif // QUANTILE_H



