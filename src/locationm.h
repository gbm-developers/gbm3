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
 CLocationM(const std::string& sType) : mtype_(sType), meps_(1e-8) {};

 CLocationM(const std::string& sType, const double& singleParam) :
  mparams_(1, singleParam), mtype_(sType), meps_(1e-8) {};

 CLocationM(const std::string& sType, const std::vector<double>& adParams) :
  mparams_(adParams), mtype_(sType), meps_(1e-8) {};

  virtual ~CLocationM() {};

  double weightedQuantile(int iN, double *adV, const double *adW, double dAlpha);

  double PsiFun(double dX);

  double LocationM(int iN, double *adX, const double *adW, double dAlpha);

private:
  std::vector<double> mparams_;
  std::string mtype_;
  double meps_;

  struct Compare
  {
    bool operator()(pair<int, double> prP, pair<int, double> prQ)
    {
      return (prP.second < prQ.second);
    }
  };
};

#endif // LOCMCGBM_H
