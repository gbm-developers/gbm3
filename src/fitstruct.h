//------------------------------------------------------------------------------
//
//  File:       fitstruct.h
//
//  Description:   struct containing parameters from fitting the gbm object
//
//  Owner:      James Hickey
//
//  History:    07/06/2016  jhickey created.
//
//------------------------------------------------------------------------------

#ifndef FITSTRUCT_H
#define FITSTRUCT_H

//------------------------------
// Includes
//------------------------------
#include "dataset.h"
#include "tree.h"

//------------------------------
// Struct definition
//------------------------------
struct FitStruct {
	//----------------------
	// Public Constructors
	//----------------------
	FitStruct():fitted_tree(),
				data_for_fit(NULL),
				training_error(0.0),
				validation_error(0.0), oobag_improvement(0.0) {
	};
	FitStruct(std::auto_ptr<CCARTTree>& tree,
			CDataset& data,
			double train_error,
			double valid_error,
			double oobag_improv):
			training_error(train_error),
			validation_error(valid_error),
			oobag_improvement(oobag_improv) {
		fitted_tree.reset(tree.release());
		data_for_fit = &data;
	};
	//----------------------
	// Public Destructors
	//----------------------
	~FitStruct() { data_for_fit = NULL; };

	//----------------------
	// Public Functions
	//----------------------
	std::auto_ptr<CCARTTree> fitted_tree;
	CDataset* data_for_fit;
	double training_error;
	double validation_error;
	double oobag_improvement;
};
#endif // FITSTRUCT_H
