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

class CLocationM {
 public:
  CLocationM(){};
  CLocationM(const std::string& kType) : mtype_(kType), meps_(1e-8){};

  CLocationM(const std::string& kType, const double& kSingleParam)
      : mparams_(1, kSingleParam), mtype_(kType), meps_(1e-8){};

  CLocationM(const std::string& kType, const std::vector<double>& kVecParams)
      : mparams_(kVecParams), mtype_(kType), meps_(1e-8){};

  virtual ~CLocationM(){};

  double WeightedQuantile(int vec_length, double* vec, const double* kWeights,
                          double alpha);

  double PsiFun(double xval);

  double LocationM(int num_data_points, double* covars, const double* kWeights,
                   double alpha);

 private:
  std::vector<double> mparams_;
  std::string mtype_;
  double meps_;

  struct Compare {
    bool operator()(pair<int, double> pair_1, pair<int, double> pair_2) {
      return (pair_1.second < pair_2.second);
    }
  };
};

#endif  // LOCMCGBM_H
