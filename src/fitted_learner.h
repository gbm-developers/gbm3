//------------------------------------------------------------------------------
//
//  File:       fitted_learner.h
//
//  Description:   class containing parameters from fitting the a tree
//
//  Owner:      James Hickey
//
//  History:    07/06/2016  jhickey created.
//
//------------------------------------------------------------------------------

#ifndef FITTEDLEARNER_H
#define FITTEDLEARNER_H

//------------------------------
// Includes
//------------------------------
#include "dataset.h"
#include "tree.h"

//------------------------------
// Class definition
//------------------------------
class FittedLearner {
public:
	//----------------------
	// Public Constructors
	//----------------------
	FittedLearner():fitted_tree_(),
				data_for_fit_(NULL),
				training_error_(0.0),
				validation_error_(0.0), oobag_improvement_(0.0) {
	};
	FittedLearner(std::auto_ptr<CCARTTree>& tree,
			CDataset& data,
			double train_error,
			double valid_error,
			double oobag_improv):
			training_error_(train_error),
			validation_error_(valid_error),
			oobag_improvement_(oobag_improv) {
		fitted_tree_.reset(tree.release());
		data_for_fit_ = &data;
	};
	//----------------------
	// Public Destructors
	//----------------------
	~FittedLearner() { data_for_fit_ = NULL; };


	//----------------------
	// Public Functions
	//----------------------
	CCARTTree* get_tree() { return fitted_tree_.get(); }
	const CCARTTree* get_tree() const { return fitted_tree_.get(); };
	CDataset* get_data_for_fit() { return data_for_fit_; } ;
	const CDataset* get_data_for_fit() const { return data_for_fit_; };
	double get_training_error() const { return training_error_; };
	double get_valid_error() const { return validation_error_; };
	double get_oobag_improv() const { return oobag_improvement_; };
	const unsigned long& get_size_of_tree() const { return fitted_tree_->size_of_tree(); };

private:
	//----------------------
	// Private Variables
	//----------------------
	std::auto_ptr<CCARTTree> fitted_tree_;
	CDataset* data_for_fit_;
	double training_error_;
	double validation_error_;
	double oobag_improvement_;
};
#endif // FITTEDLEARNER_H
