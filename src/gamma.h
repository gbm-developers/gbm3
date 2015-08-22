#ifndef GAMMA_H
#define GAMMA_H

#include <Rmath.h>
#include "distribution.h"

class CGamma : public CDistribution
{

public:

    CGamma();

    virtual ~CGamma();

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
    vector<double> vecdMax;
    vector<double> vecdMin;
};

#endif // GAMMA_H



