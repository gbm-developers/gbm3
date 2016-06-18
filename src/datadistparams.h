//------------------------------------------------------------------------------
//
//  File:       datadistparams.h
//
//  Description: struct used to initialize data/dist container
//
//  Owner:      James Hickey
//
//  History:    06/06/2016  James Hickey created.
//
//------------------------------------------------------------------------------

#ifndef DATADISTPARAMS_H
#define DATADISTPARAMS_H

//------------------------------
// Includes
//------------------------------
#include "gbm_exception.h"
#include "parallel_details.h"

#include <Rcpp.h>

//------------------------------
// class definition
//------------------------------
class DataDistParams {
 public:
  //----------------------
  // Public Constructors
  //----------------------
  //-----------------------------------
  // Function: DataDistParams
  //
  // Returns: None
  //
  // Description: Constructor for DataDistParams struc.
  //
  // Parameters:
  //  response  - SEXP containing the response of each data-point - accessed via
  //				double ptr
  //  offset_vec - SEXP containing the offset applied to each response -
  //  accessed via
  //				double ptr
  //  covariates - SEXP containing the predictor values - becomes
  //  Rcpp::NumericMatrix
  //  covar_order - SEXP containing the order of predictor values to
  //  			be used in GBM formula - accessed via int ptr.
  //  sorted_vec - SEXP indicating the ordering of observations for CoxPH
  //				 , this is stored an int ptr in CoxPH.
  //  strata_vec - SEXP indicating which strata observations are in -
  //				 this is used with CoxPH as an int ptr.
  //  obs_weight  - SEXP containing weights to be used in fitting
  //  				process - accessed via double ptr.
  //  misc - SEXP list object containing distribution dependent data
  //		   , depending on the distribution this may be a double or int.
  // prior_coeff_var - SEXP containing a double specifying a prior node
  //					prediction value for the CoxPH model.
  //  row_to_obs_id - SEXP storing integer ids mapping each row to
  //          		its observation - becomes Rcpp::IntegerVector.
  //  var_classes  - R object (SEXP) containing the variable classes,
  //				  this is stored as a Rcpp::IntegerVector in
  // dataset.cpp.
  //  monotonicity_vec - SEXP containing +-1/0 indicating whether the
  //  					variables are:
  //  monotone increasing (+1), decreasing (-1) or arbitrary (0) with the
  //  response
  // 					 variable. - accessed via int ptr
  //  dist_family - SEXP specifying distribution to instantiate - string.
  //  fraction_in_bag - SEXP specifying the fraction of observations to be
  //  bagged -
  //  				   which is a double.
  //  num_rows_in_training - SEXP containing ulong specify the number of data
  //  points in
  //							training set
  //  unique_training_obs - SEXP containing ulong specifying number of unique
  //  observations
  //						  that make up training set
  //  number_offeatures - SEXP containing ulong specifying number of features
  //  for use in
  //						tree growing.
  //  parallel - parallelization-related constants
  //-----------------------------------
  DataDistParams(SEXP response, SEXP offset_vec, SEXP covariates,
                 SEXP covar_order, SEXP sorted_vec, SEXP strata_vec,
                 SEXP obs_weight, SEXP misc, SEXP prior_coeff_var,
                 SEXP row_to_obs_id, SEXP var_classes, SEXP monotonicity_vec,
                 SEXP dist_family, SEXP fraction_inbag,
                 SEXP num_rows_in_training, SEXP unique_training_obs,
                 SEXP number_offeatures, const parallel_details& parallel)
      : response(response),
        observationids(row_to_obs_id),
        misc(misc),
        parallel(parallel) {
    sorted = sorted_vec;
    strata = strata_vec;
    offset = offset_vec;
    xvalues = covariates;
    xorder = covar_order;
    variable_weight = obs_weight;
    variable_num_classes = var_classes;
    variable_monotonicity = monotonicity_vec;
    num_trainrows = Rcpp::as<unsigned long>(num_rows_in_training);
    num_trainobservations = Rcpp::as<unsigned long>(unique_training_obs);
    num_features = Rcpp::as<unsigned long>(number_offeatures);
    bagfraction = Rcpp::as<double>(fraction_inbag);
    prior_coefficient_variation = Rcpp::as<double>(prior_coeff_var);

    // Set up distribution family
    family = Rcpp::as<std::string>(dist_family);
    if (family.empty()) {
      throw gbm_exception::Failure(
          "configStructs - Can't specify IR metric as family not initialized.");
    }
  }

  //-------------------
  // Public Variables
  //-------------------
  Rcpp::NumericMatrix response;
  Rcpp::IntegerVector observationids;
  Rcpp::List misc;
  parallel_details parallel;
  SEXP sorted;
  SEXP strata;
  SEXP offset;
  SEXP xvalues;
  SEXP xorder;
  SEXP variable_weight;
  SEXP variable_num_classes;
  SEXP variable_monotonicity;
  unsigned long num_trainrows;
  unsigned long num_trainobservations;
  unsigned long num_features;
  double bagfraction;
  double prior_coefficient_variation;
  std::string family;
};
#endif  // DATADISTPARAMS_H
