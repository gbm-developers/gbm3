//------------------------------------------------------------------------------
//  GBM alteration by Daniel Edwards
//  File:       locationm.h
//
//  History:    27/3/2008 created
//
//------------------------------------------------------------------------------

#ifndef LOCMCGBM_H
#define LOCMCGBM_H

#include <string>
#include <utility>
#include <vector>
#include <algorithm>
#include <R.h>

using namespace std;

class CLocationM
{

public:

 CLocationM(){};
 CLocationM(const std::string& sType) : msType(sType), mdEps(1e-8) {};

 CLocationM(const std::string& sType, const double& singleParam) :
  madParams(1, singleParam), msType(sType), mdEps(1e-8) {};

 CLocationM(const std::string& sType, const std::vector<double>& adParams) :
  madParams(adParams), msType(sType), mdEps(1e-8) {};

  virtual ~CLocationM() {};

  double weightedQuantile(int iN, double *adV, const double *adW, double dAlpha);

  double PsiFun(double dX);

  double LocationM(int iN, double *adX, const double *adW, double dAlpha);

private:
  std::vector<double> madParams;
  std::string msType;
  double mdEps;

  struct Compare
  {
    bool operator()(pair<int, double> prP, pair<int, double> prQ)
    {
      return (prP.second < prQ.second);
    }
  };
};

#endif // LOCMCGBM_H
