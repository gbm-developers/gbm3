//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//  File:       poisson.h
//
//  License:    GNU GPL (version 2 or later)
//
//  Contents:   poisson object
//
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef POISSON_H
#define POISSON_H

#include <Rmath.h>
#include "distribution.h"

class CPoisson : public CDistribution
{

public:

    CPoisson();

    virtual ~CPoisson();

    
    void ComputeWorkingResponse(const double *adY,
				const double *adMisc,
				const double *adOffset,
				const double *adF,
				double *adZ,
				const double *adWeight,
				const bag& afInBag,
				unsigned long nTrain);

    double Deviance(const double *adY,
                    const double *adMisc,
                    const double *adOffset,
                    const double *adWeight,
                    const double *adF,
                    unsigned long cLength);

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
    vector<double> vecdNum;
    vector<double> vecdDen;
    vector<double> vecdMax;
    vector<double> vecdMin;
};

#endif // POISSON_H



