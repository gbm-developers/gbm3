//------------------------------------------------------------------------------
//  GBM alteration by Daniel Edwards
//  File:       locationm.h
//
//  History:    27/3/2008 created
//        
//------------------------------------------------------------------------------

#ifndef LOCMCGBM_H
#define LOCMCGBM_H

#include <utility>
#include <vector>
#include <algorithm>
#include <R.h>

using namespace std;

class CLocationM
{

public:

    CLocationM(const char *sType, int iN, double *adParams);

    virtual ~CLocationM();

	double Median(int iN, double *adV, double *adW);

	double PsiFun(double dX);

	double LocationM(int iN, double *adX, double *adW);

private:
	double *madParams;
	const char *msType;
	double mdEps;

    struct comp{

	    bool operator()(pair<int, double> prP, pair<int, double> prQ)
		{
		    return (prP.second < prQ.second);
		}
	};
};

#endif // LOCMCGBM_H



