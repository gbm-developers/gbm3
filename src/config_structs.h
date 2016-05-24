//------------------------------------------------------------------------------
//
//  File:       configStructs.h
//
//  Description:   Structs used to wrap up initialization parameters for
//				both the weak learner container and data/dist container
//
//  Owner:      James Hickey
//
//  History:    31/03/2016  James Hickey created.
//
//------------------------------------------------------------------------------

#ifndef CONFIGSTRUCTS_H
#define CONFIGSTRUCTS_H

//------------------------------
// Includes
//------------------------------
#include "gbmexcept.h"
#include <Rcpp.h>

//------------------------------
// Struct definitions
//------------------------------
class DataDistParams
{
public:
	DataDistParams(SEXP response,
		    SEXP offset_vec,
		    SEXP covariates,
		    SEXP covar_order,
		    SEXP sorted_vec,
		    SEXP strata_vec,
		    SEXP obs_weight,
		    SEXP misc,
		    SEXP prior_coeff_var,
		    SEXP row_to_obs_id,
		    SEXP var_classes,
		    SEXP monotonicity_vec,
		    SEXP dist_family,
		    SEXP fraction_inbag,
		    SEXP num_rows_in_training,
		    SEXP unique_training_obs,
		    SEXP number_offeatures):
		    response(response),
		    patientids(row_to_obs_id),
		    misc(misc)
{
		sorted = sorted_vec;
		strata = strata_vec;
		offset = offset_vec;
		xvalues = covariates;
		xorder = covar_order;
		variable_weight =	obs_weight;
		variable_num_classes = var_classes;
		variable_monotonicity = monotonicity_vec;
		num_trainrows = Rcpp::as<int>(num_rows_in_training);
		num_trainobservations = Rcpp::as<int>(unique_training_obs);
		num_features = Rcpp::as<int>(number_offeatures);
		bagfraction = Rcpp::as<double>(fraction_inbag);
		prior_coefficient_variation = Rcpp::as<double>(prior_coeff_var);
		family = Rcpp::as<std::string>(dist_family);
		irmeasure = NULL;

}
	Rcpp::NumericMatrix response;
	Rcpp::IntegerVector patientids;
	Rcpp::List misc;
	SEXP sorted;
	SEXP strata;
	SEXP offset;
	SEXP xvalues;
	SEXP xorder;
	SEXP variable_weight;
	SEXP variable_num_classes;
	SEXP variable_monotonicity;
	int num_trainrows;
	int num_trainobservations;
	int num_features;
	double bagfraction;
	double prior_coefficient_variation;
	std::string family;
	const char* irmeasure;
};

struct TreeParams
{
  int depth;
  int min_obs_in_node;
  double shrinkage;
  int num_trainrows;
  int numberdatacolumns;
};

//------------------------------
// Class definition
//------------------------------
class ConfigStructs
{
public:
	//----------------------
	// Public Constructors
	//----------------------
	ConfigStructs(SEXP response,
		    SEXP offset_vec,
		    SEXP covariates,
		    SEXP covar_order,
		    SEXP sorted_vec,
		    SEXP strata_vec,
		    SEXP obs_weight,
		    SEXP misc,
		    SEXP prior_coeff_var,
		    SEXP row_to_obs_id,
		    SEXP var_classes,
		    SEXP monotonicity_vec,
		    SEXP dist_family,
		    SEXP num_trees,
		    SEXP tree_depth,
		    SEXP min_num_node_obs,
		    SEXP shrinkageconstant,
		    SEXP fraction_inbag,
		    SEXP num_rows_in_training,
		    SEXP unique_training_obs,
		    SEXP number_offeatures):

		    dataconfig_(response,
				    offset_vec,
				    covariates,
				    covar_order,
				    sorted_vec,
				    strata_vec,
				    obs_weight,
				    misc,
				    prior_coeff_var,
				    row_to_obs_id,
				    var_classes,
				    monotonicity_vec,
				    dist_family,
				    fraction_inbag,
				    num_rows_in_training,
				    unique_training_obs,
				    number_offeatures)
	{
		// Initialize DataDistParams
		InitIRMeasure();

		// Initialize Tree Parameters
		treeconfig_.depth = Rcpp::as<int>(tree_depth);
		treeconfig_.min_obs_in_node = Rcpp::as<int>(min_num_node_obs);
		treeconfig_.shrinkage = Rcpp::as<double>(shrinkageconstant);
		treeconfig_.num_trainrows = Rcpp::as<int>(num_rows_in_training);

		Rcpp::NumericMatrix covar_mat(covariates);
		treeconfig_.numberdatacolumns =  covar_mat.ncol();
	};


	//---------------------
	// Public destructor
	//---------------------
	~ConfigStructs(){};

	//---------------------
	// Public Methods
	//---------------------
	DataDistParams& get_data_config() { return dataconfig_; };
	TreeParams& get_tree_config() { return treeconfig_; };

private:
	//-------------------
	// Private Variables
	//-------------------
	DataDistParams dataconfig_;
	TreeParams treeconfig_;

	//-------------------
	// Private Methods
	//-------------------
	void InitIRMeasure()
	{
		// Check family specified
		if(dataconfig_.family.empty())
		{
			throw GBM::Failure("configStructs - Can't specify IR metric as family not initialized.");
		}
		// Get szIRMeasure and throw appropriate errors
		if(0 == dataconfig_.family.compare(0, 8, "pairwise"))
		{
			std::size_t offset_tomeasure = dataconfig_.family.find("_");
			if(offset_tomeasure == std::string::npos)
			{
				throw GBM::Failure("Unable to locate IR metric required for pairwise");
			}
			const char* kIrMeasure = dataconfig_.family.c_str() + offset_tomeasure + 1;
			dataconfig_.irmeasure = kIrMeasure;
			dataconfig_.family = "pairwise";
		}
		else
		{
			dataconfig_.irmeasure = "";
		}


	}

};
#endif // CONFIGSTRUCTS_H
