//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//
//  File:       bernoulli.h
//
//  License:    GNU GPL (version 2 or later)
//
//  Contents:   bernoulli object
//
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef BERNOULLI_H
#define BERNOULLI_H

#include "distribution.h"
#include "buildinfo.h"

class CBernoulli : public CDistribution
{

public:

    CBernoulli();

    virtual ~CBernoulli();

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
    bool fCappedPred;
};

#endif // BERNOULLI_H



