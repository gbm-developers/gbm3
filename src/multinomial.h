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
#include "distribution.h"
#include "locationm.h"

class CMultinomial : public CDistribution
{

public:

    CMultinomial(int cNumClasses, int cRows);

    virtual ~CMultinomial();
    GBMRESULT UpdateParams(double *adF,
                           double *adOffset,
                           double *adWeight,
                           unsigned long cLength);

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
   unsigned long mcNumClasses;
   unsigned long mcRows;
   double *madProb; 
};

#endif // KMULTICGBM_H


