//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//
//  File:       dataset.h
//
//  License:    GNU GPL (version 2 or later)
//
//  Contents:   Dataset class
//
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef DATASET_H
#define DATASET_H

#include <vector>

#include <Rcpp.h>

typedef std::vector<int> bag;

#include "buildinfo.h"
#include "gbmexcept.h"

namespace {
  inline bool has_value(const Rcpp::NumericVector& x) {
    return !( (x.size() == 1) && (ISNA(x[0])));
  }

  std::ptrdiff_t shuffler(std::ptrdiff_t n) {
    return n * unif_rand();
  }
}
  

class CDataset
{
public:
  CDataset(SEXP radY, SEXP radOffset, SEXP radX, SEXP raiXOrder,
           SEXP radWeight, SEXP radMisc,
           SEXP racVarClasses, SEXP ralMonotoneVar) :
  adY(radY), adOffset(radOffset), adWeight(radWeight), adMisc(radMisc),
    adX(radX),
    acVarClasses(racVarClasses), alMonotoneVar(ralMonotoneVar),
    aiXOrder(raiXOrder),
    fHasMisc(has_value(adMisc)), fHasOffset(has_value(adOffset)) {

    if (adX.ncol() != alMonotoneVar.size()) {
      throw GBM::invalid_argument("shape mismatch (monotone does not match data)");
    }
    
    if (adX.ncol() != acVarClasses.size()) {
      throw GBM::invalid_argument("shape mismatch (var classes does not match daa)");
    }

  };

  virtual ~CDataset()  {};

  typedef std::vector<int> index_vector;
  
  int nrow() const {
    return adX.nrow();
  }

  int ncol() const {
    return adX.ncol();
  }

  double* y_ptr() {
    return adY.begin();
  }

  const double* y_ptr() const {
    return adY.begin();
  }

  
  double* offset_ptr(bool require=true) {
    if (has_offset()) {
      return adOffset.begin();
    } else {
      if (require) {
        throw GBM::failure("You require a genuine offset, and don't have one.");
      } else {
        return 0;
      }
    }
  }    

  const double* offset_ptr(bool require=true) const {
    if (has_offset()) {
      return adOffset.begin();
    } else {
      if (require) {
        throw GBM::failure("You require a genuine offset, and don't have one.");
      } else {
        return 0;
      }
    }
  }    

  double* weight_ptr() {
    return adWeight.begin();
  }

  const double* weight_ptr() const {
    return adWeight.begin();
  }

  double* misc_ptr(bool require=true) {
    if (has_misc()) {
      return adMisc.begin();
    } else {
      if (require) {
        throw GBM::failure("You require genuine misc, and don't have it.");
      } else {
        return 0;
      }
    }
  }

  const double* misc_ptr(bool require=true) const {
    if (has_misc()) {
      return adMisc.begin();
    } else {
      if (require) {
        throw GBM::failure("You require genuine misc, and don't have it.");
      } else {
        return 0;
      }
    }
  }
  
  double* x_ptr() {
    return adX.begin();
  }

  const double* x_ptr() const {
    return adX.begin();
  }

  int varclass(int ind) const {
    return acVarClasses[ind];
  }

  int* monotone_ptr() {
    return alMonotoneVar.begin();
  }

  const int* monotone_ptr() const {
    return alMonotoneVar.begin();
  }
  
  int* order_ptr() {
    return aiXOrder.begin();
  }

  const int* order_ptr() const {
    return aiXOrder.begin();
  }

  bool has_misc() const {
    return fHasMisc;
  }

  bool has_offset() const {
    return fHasOffset;
  }

  double x_value(const int row, const int col) const {
    return adX(row, col);
  }

  index_vector random_order() const {
    index_vector result(ncol());
    // fill the vector
    for (index_vector::size_type ind=0; ind!=result.size(); ++ind) {
      result[ind] = ind;
    }
    // and now shuffle
    std::random_shuffle(result.begin(), result.end(), shuffler);
    // and return
    return result;
  }
  
 private:
    
  Rcpp::NumericVector adY, adOffset, adWeight, adMisc;
  Rcpp::NumericMatrix adX;
  Rcpp::IntegerVector acVarClasses, alMonotoneVar, aiXOrder;

  bool fHasMisc;
  bool fHasOffset;
  
};

#endif // DATASET_H


