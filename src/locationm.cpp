//------------------------------------------------------------------------------
//  GBM alteration by Daniel Edwards
//  File:       locationm.cpp
//
//  Purpose:    Class to provide methods to calculate the location M-estimates
//              of a variety of functions
//
//  History:    31/03/2008 created
//
//------------------------------------------------------------------------------

#include "locationm.h"
#include "gbm_exception.h"

#include <algorithm>
#include <Rcpp.h>

/////////////////////////////////////////////////
// weightedQuantile
//
// Function to return the weighted quantile of
// a vector of a given length
//
// Parameters: iN     - Length of vector
//             adV    - Vector of doubles
//             adW    - Array of weights
//             dAlpha - Quantile to calculate (0.5 for median)
//
// Returns :   Weighted quantile
/////////////////////////////////////////////////
double CLocationM::WeightedQuantile(int vec_length, double* vec,
                                    const double* kWeights, double alpha) {
  // Local variables
  int ii, med_idx;
  std::vector<double> vec_w;
  std::vector<pair<int, double> > vec_v;
  double cum_sum, wsum, med;

  // Check the vector size
  if (vec_length == 0) {
    return 0.0;
  } else if (vec_length == 1) {
    return vec[0];
  }

  // Create vectors containing the values and weights
  vec_v.resize(vec_length);
  for (ii = 0; ii < vec_length; ii++) {
    vec_v[ii] = make_pair(ii, vec[ii]);
  }

  // Sort the vector
  std::stable_sort(vec_v.begin(), vec_v.end(), Compare());

  // Sort the weights correspondingly and calculate their sum
  vec_w.resize(vec_length);
  wsum = 0.0;
  for (ii = 0; ii < vec_length; ii++) {
    vec_w[ii] = kWeights[vec_v[ii].first];
    wsum += kWeights[ii];
  }

  // Get the first index where the cumulative weight is >=0.5
  med_idx = -1;
  cum_sum = 0.0;
  while (cum_sum < alpha * wsum) {
    med_idx++;
    cum_sum += vec_w[med_idx];
  }

  // Get the index of the next non-zero weight
  int iNextNonZero = vec_length;
  for (ii = (vec_length - 1); ii > med_idx; ii--) {
    if (vec_w[ii] > 0) {
      iNextNonZero = ii;
    }
  }

  // Use this index unless the cumulative sum is exactly alpha
  if (iNextNonZero == vec_length || cum_sum > alpha * wsum) {
    med = vec_v[med_idx].second;
  } else {
    med = alpha * (vec_v[med_idx].second + vec_v[iNextNonZero].second);
  }

  return med;
}

/////////////////////////////////////////////////
// PsiFun
//
// Function to calculate the psi of the supplied
// value, given the type of function to use and
// the supplied parameters
//
// Parameters: dX - Value
//
// Returns :   Psi(X)
/////////////////////////////////////////////////
double CLocationM::PsiFun(double xval) {
  /*
    why we are checking this here rather than at construction
    is entirely unknown...
   */
  if (mtype_ == "tdist") {
    return xval / (mparams_[0] + (xval * xval));
  }

  throw gbm_exception::Failure("Function type " + mtype_ + "not known.");
}

/////////////////////////////////////////////////
// LocationM
//
// Function to calculate location M estimate for
// the supplied weighted data, with the psi-function
// type and parameters specified in this class
//
// Parameters: iN  - Number of data points
//             adX - Data vector
//             adW - Weight vector
//             dAlpha - Quantile to calculate (0.5 for median)
//
// Returns :   Location M-Estimate of (X, W)
/////////////////////////////////////////////////
double CLocationM::LocationM(int num_data_points, double* covars,
                             const double* kWeights, double alpha) {
  // Local variables
  int ii;

  // Get the initial estimate of location
  double beta0 = WeightedQuantile(num_data_points, covars, kWeights, alpha);

  // Get the initial estimate of scale
  std::vector<double> diff_vec(num_data_points);
  for (ii = 0; ii < num_data_points; ii++) {
    diff_vec[ii] = fabs(covars[ii] - beta0);
  }

  double scale0 =
      1.4826 * WeightedQuantile(num_data_points, &diff_vec[0], kWeights, alpha);
  scale0 = R::fmax2(scale0, meps_);

  // Loop over until the error is low enough
  double error = 1.0;
  int counter = 0;

  while (counter < 50) {
    double sum_w_x = 0.0;
    double sum_w = 0.0;
    for (ii = 0; ii < num_data_points; ii++) {
      double dt = fabs(covars[ii] - beta0) / scale0;
      dt = R::fmax2(dt, meps_);
      double dwt = kWeights[ii] * PsiFun(dt) / dt;

      sum_w_x += dwt * covars[ii];
      sum_w += dwt;
    }

    double beta = beta0;
    if (sum_w > 0) {
      beta = sum_w_x / sum_w;
    }

    error = fabs(beta - beta0);
    if (error > meps_) {
      error /= fabs(beta0);
    }
    beta0 = beta;

    if (error < meps_) {
      counter = 100;
    } else {
      counter++;
    }
  }

  return beta0;
}
