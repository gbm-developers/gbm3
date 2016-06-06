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

#ifndef TREEPARAMS_H
#define TREEPARAMS_H

//------------------------------
// Includes
//------------------------------
#include "gbm_exception.h"
#include <Rcpp.h>

//------------------------------
// Struct definitions
//------------------------------
struct TreeParams {

  //----------------------
  // Public Constructors
  //----------------------
  TreeParams(SEXP tree_depth, SEXP min_num_node_obs,
		  	 SEXP shrinkageconstant, SEXP num_rows_in_training) {

	  depth = Rcpp::as<unsigned long>(tree_depth);
	  min_obs_in_node = Rcpp::as<unsigned long>(min_num_node_obs);
	  shrinkage = Rcpp::as<double>(shrinkageconstant);
	  num_trainrows = Rcpp::as<unsigned long>(num_rows_in_training);
  }

  //----------------------
  // Public Constructors
  //----------------------
  unsigned long depth;
  unsigned long min_obs_in_node;
  double shrinkage;
  unsigned long num_trainrows;
};
#endif // TREEPARAMS_H
