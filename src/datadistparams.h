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
#include <Rcpp.h>

//------------------------------
// class definition
//------------------------------
class DataDistParams {
public:
  //----------------------
  // Public Constructors
  //----------------------
  DataDistParams(SEXP response, SEXP offset_vec, SEXP covariates,
                 SEXP covar_order, SEXP sorted_vec, SEXP strata_vec,
                 SEXP obs_weight, SEXP misc, SEXP prior_coeff_var,
                 SEXP row_to_obs_id, SEXP var_classes, SEXP monotonicity_vec,
                 SEXP dist_family, SEXP fraction_inbag,
                 SEXP num_rows_in_training, SEXP unique_training_obs,
                 SEXP number_offeatures)
      : response(response), observationids(row_to_obs_id), misc(misc) {
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
    family = Rcpp::as<std::string>(dist_family);

    // Set up IR measure
    irmeasure = NULL;
    InitIRMeasure();

  }

  //-------------------
  // Public Variables
  //-------------------
  Rcpp::NumericMatrix response;
  Rcpp::IntegerVector observationids;
  Rcpp::List misc;
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
  const char* irmeasure;

private:
  //-------------------
  // Private Methods
  //-------------------
  void InitIRMeasure();
};
#endif  // DATADISTPARAMS_H
