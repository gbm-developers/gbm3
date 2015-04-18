#ifndef TWEEDIE_H
#define TWEEDIE_H

#include <Rmath.h>
#include "distribution.h"

class CTweedie : public CDistribution
{

public:

    CTweedie(double dPower) : dPower(dPower) {};

    virtual ~CTweedie();

    void ComputeWorkingResponse(const double *adY,
				const double *adMisc,
				const double *adOffset,
				const double *adF,
				double *adZ,
				const double *adWeight,
				const bag& afInBag,
				unsigned long nTrain,
				int cIdxOff);

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
			 const double *adFadj,
			 int cIdxOff);

    double Deviance(const double *adY,
                    const double *adMisc,
                    const double *adOffset,
                    const double *adWeight,
                    const double *adF,
                    unsigned long cLength,
		    int cIdxOff);
    
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
    double dPower;
};

#endif // TWEEDIE_H

