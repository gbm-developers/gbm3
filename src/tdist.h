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

    CTDist(double adNu);

    virtual ~CTDist();

	GBMRESULT UpdateParams(double *adF,
	                       double *adOffset,
						   double *adWeight,
	                       unsigned long cLength)
	{ 
		return GBM_OK;
	};

    GBMRESULT ComputeWorkingResponse(double *adY,
                                     double *adMisc,
                                     double *adOffset,
                                     double *adF, 
                                     double *adZ, 
                                     double *adWeight,
                                     bool *afInBag,
                                     unsigned long nTrain,
	                                 int cIdxOff);
 
    GBMRESULT InitF(double *adY, 
                    double *adMisc,
                    double *adOffset,
                    double *adWeight,
                    double &dInitF, 
                    unsigned long cLength);

    GBMRESULT FitBestConstant(double *adY,
                              double *adMisc,
                              double *adOffset,
                              double *adW,
                              double *adF,
                              double *adZ,
                              unsigned long *aiNodeAssign,
                              unsigned long nTrain,
                              VEC_P_NODETERMINAL vecpTermNodes,
                              unsigned long cTermNodes,
                              unsigned long cMinObsInNode,
                              bool *afInBag,
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
                          bool *afInBag,
                          double dStepSize,
                          unsigned long nTrain);

private:
	double mdNu;
	CLocationM *mpLocM;
};

#endif // TDISTCGBM_H



