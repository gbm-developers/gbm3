//------------------------------------------------------------------------------
//
//  File:       gbm_engine.h
//
//  Description:   Header file for Gradient Boosting Engine.
//
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef GBMENGINE_H
#define GBMENGINE_H

//------------------------------
// Includes
//------------------------------
#include "config_structs.h"
#include "gbm_datacontainer.h"
#include "tree.h"
#include <memory>
#include <Rcpp.h>
#include <vector>

//------------------------------
// Class definition
//------------------------------
class CGBM
{
public:
	//----------------------
	// Public Constructors
	//----------------------
    CGBM(ConfigStructs& gbmparams);

	//---------------------
	// Public destructor
	//---------------------
    ~CGBM();

	//---------------------
	// Public Functions
	//---------------------
    void FitLearner(double* func_estimate,
		 double &training_error,
		 double &validation_error,
		 double &outofbag_improv);

    void GbmTransferTreeToRList(int* splitvar,
			     double* splitvalues,
			     int* leftnodes,
			     int* rightnodes,
			     int* missingnodes,
			     double* error_reduction,
			     double* weights,
			     double* predictions,
			     VecOfVectorCategories &splitcodes_vec,
			     int prev_categorical_splits);

    const long size_of_fitted_tree() const{ return tree_.size_of_tree(); }
    double initial_function_estimate() { return datacontainer_.InitialFunctionEstimate(); };

private:
	//-------------------
	// Private Variables
	//-------------------
    CGBMDataContainer datacontainer_;
    CCARTTree tree_;
    
    // Residuals and adjustments to function estimate
    std::vector<double> residuals_;

};

#endif // GBMENGINE_H



