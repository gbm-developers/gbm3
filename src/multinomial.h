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
                                     int *afInBag,
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
   unsigned long mcNumClasses;
   unsigned long mcRows;
   std::vector<double> madProb;
};

#endif // KMULTICGBM_H


