//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//
//  File:       adaboost.h
//
//  License:    GNU GPL (version 2 or later)
//
//  Contents:   Object for fitting for the AdaBoost loss function
//
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef ADABOOST_H
#define ADABOOST_H

#include "distribution.h"

class CAdaBoost : public CDistribution
{

public:

    CAdaBoost();

    virtual ~CAdaBoost();

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
    vector<double> vecdNum;
    vector<double> vecdDen;
};

#endif // ADABOOST_H


