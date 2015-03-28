#ifndef GAMMA_H
#define GAMMA_H

#include <Rmath.h>
#include "distribution.h"

class CGamma : public CDistribution
{

public:

    CGamma();

    virtual ~CGamma();

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
    vector<double> vecdMax;
    vector<double> vecdMin;
};

#endif // GAMMA_H



