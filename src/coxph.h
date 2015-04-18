//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//
//  File:       coxph.h
//
//  License:    GNU GPL (version 2 or later)
//
//  Contents:   Cox proportional hazard object
//
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------


#ifndef COXPH_H
#define COXPH_H

#include "distribution.h"
#include "matrix.h"

class CCoxPH : public CDistribution
{

public:

    CCoxPH();

    virtual ~CCoxPH();

    void ComputeWorkingResponse(const double *adT,
				const double *adDelta,
				const double *adOffset,
				const double *adF,
				double *adZ,
				const double *adWeight,
				const bag& afInBag,
				unsigned long nTrain,
				int cIdxOff);

    void InitF(const double *adT,
	       const double *adDelta,
	       const double *adOffset,
	       const double *adWeight,
	       double &dInitF,
	       unsigned long cLength);
    
    void FitBestConstant(const double *adT,
			 const double *adDelta,
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
    
    double Deviance(const double *adT,
                    const double *adDelta,
                    const double *adOffset,
                    const double *adWeight,
                    const double *adF,
                    unsigned long cLength,
		    int cIdxOff);

    double BagImprovement(const double *adT,
                          const double *adDelta,
                          const double *adOffset,
                          const double *adWeight,
                          const double *adF,
                          const double *adFadj,
                          const bag& afInBag,
                          double dStepSize,
                          unsigned long nTrain);


private:
    vector<double> vecdP;
    vector<double> vecdRiskTot;
    vector<double> vecdG;
    vector<unsigned long> veciK2Node;
    vector<unsigned long> veciNode2K;

    matrix<double> matH;
    matrix<double> matHinv;
};

#endif // COXPH_H



