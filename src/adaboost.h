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

    void ComputeWorkingResponse(double *adY,
				double *adMisc,
				double *adOffset,
				double *adWeight,
				double *adF,
				double *adZ,
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
			 const std::vector<unsigned long>& aiNodeAssign,
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
    vector<double> vecdNum;
    vector<double> vecdDen;
};

#endif // ADABOOST_H


