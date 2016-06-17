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
#include "parallel_details.h"
#include <Rcpp.h>

//------------------------------
// Struct definitions
//------------------------------
struct TreeParams {
  //----------------------
  // Public Constructors
  //----------------------
  //-----------------------------------
  // Function: TreeParams
  //
  // Returns: None
  //
  // Description: Constructor for TreeParams struc.
  //
  // Parameters:
  //  tree_depth - SEXP specifying the maximum depth of each tree - ulong.
  //  min_num_node_obs - SEXP specifying the minimum number of obs.
  //					a node must have - ulong.
  //  shrinkageconstant - SEXP defining the shrinkage applied to
  //				  each tree fit - double.
  //  num_rows_in_training - SEXP containing ulong specify the number of data
  //  points in
  //							training set
  //  parallel - parallelization-related constants
  //-----------------------------------

  TreeParams(SEXP tree_depth, SEXP min_num_node_obs, SEXP shrinkageconstant,
             SEXP num_rows_in_training, const parallel_details& parallel)
      : depth(Rcpp::as<unsigned long>(tree_depth)),
        min_obs_in_node(Rcpp::as<unsigned long>(min_num_node_obs)),
        shrinkage(Rcpp::as<double>(shrinkageconstant)),
        num_trainrows(Rcpp::as<unsigned long>(num_rows_in_training)),
        parallel(parallel) {}

  //----------------------
  // Public Constructors
  //----------------------
  unsigned long depth;
  unsigned long min_obs_in_node;
  double shrinkage;
  unsigned long num_trainrows;
  parallel_details parallel;
};
#endif  // TREEPARAMS_H
