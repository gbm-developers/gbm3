//-----------------------------------
//
// File: dataset.cpp
//
// Description: implementation of datadistparams
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "dataset.h"
#include <Rcpp.h>

void DataDistParams::InitIRMeasure() {
	// Check family specified
	if (family.empty()) {
	  throw gbm_exception::Failure(
		  "configStructs - Can't specify IR metric as family not initialized.");
	}
	// Get szIRMeasure and throw appropriate errors
	if (0 == family.compare(0, 8, "pairwise")) {
	  std::size_t offset_tomeasure = family.find("_");
	  if (offset_tomeasure == std::string::npos) {
		throw gbm_exception::Failure(
			"Unable to locate IR metric required for pairwise");
	  }
	  const char* kIrMeasure =
		  family.c_str() + offset_tomeasure + 1;
	  irmeasure = kIrMeasure;
	  family = "pairwise";
	} else {
	  irmeasure = "";
	}
}
