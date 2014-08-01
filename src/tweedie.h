#ifndef TWEEDIE_H
#define TWEEDIE_H

#include <Rmath.h>
#include "distribution.h"

class Ctweedie : public CDistribution
{

public:

    Ctweedie();

    virtual ~Ctweedie();

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
                                     double *adWeight,
                                     double *adF,
                                     double *adZ,
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
    vector<double> vecdNum;
    vector<double> vecdDen;
	vector<double> vecdMax;
    vector<double> vecdMin;
};

#endif // TWEEDIE_H

