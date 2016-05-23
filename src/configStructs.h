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
	DataDistParams(SEXP radY,
		    SEXP radOffset,
		    SEXP radX,
		    SEXP raiXOrder,
		    SEXP rSorted,
		    SEXP rStrata,
		    SEXP radWeight,
		    SEXP radMisc,
		    SEXP rPriorCoeff,
		    SEXP rPatientId,
		    SEXP racVarClasses,
		    SEXP ralMonotoneVar,
		    SEXP rszFamily,
		    SEXP rdBagFraction,
		    SEXP rcTrain,
		    SEXP rcTrainPatients,
		    SEXP rcFeatures):
		    respY(radY),
		    patId(rPatientId),
		    misc(radMisc)
{
		sorted = rSorted;
		strata = rStrata;
		offset = radOffset;
		xValues = radX;
		xOrder = raiXOrder;
		varWeight =	radWeight;
		varClasses = racVarClasses;
		monotoneVar = ralMonotoneVar;
		cTrain = Rcpp::as<int>(rcTrain);
		cTrainPatients = Rcpp::as<int>(rcTrainPatients);
		cFeatures = Rcpp::as<int>(rcFeatures);
		dBagFraction = Rcpp::as<double>(rdBagFraction);
		priorCoeffVar = Rcpp::as<double>(rPriorCoeff);
		family = Rcpp::as<std::string>(rszFamily);
		szIRMeasure = NULL;

}
	Rcpp::NumericMatrix respY;
	Rcpp::IntegerVector patId;
	Rcpp::List misc;
	SEXP sorted;
	SEXP strata;
	SEXP offset;
	SEXP xValues;
	SEXP xOrder;
	SEXP varWeight;
	SEXP varClasses;
	SEXP monotoneVar;
	int cTrain;
	int cTrainPatients;
	int cFeatures;
	double dBagFraction;
	double priorCoeffVar;
	std::string family;
	const char* szIRMeasure;
};

struct TreeParams
{
  int cDepth;
  int cMinObsInNode;
  double dShrinkage;
  int cTrain;
  int numColData;
};

//------------------------------
// Class definition
//------------------------------
class configStructs
{
public:
	//----------------------
	// Public Constructors
	//----------------------
	configStructs(SEXP radY,
		    SEXP radOffset,
		    SEXP radX,
		    SEXP raiXOrder,
		    SEXP rSorted,
		    SEXP rStrata,
		    SEXP radWeight,
		    SEXP radMisc,
		    SEXP rPriorCoeff,
		    SEXP rPatientId,
		    SEXP racVarClasses,
		    SEXP ralMonotoneVar,
		    SEXP rszFamily,
		    SEXP rcTrees,
		    SEXP rcDepth,
		    SEXP rcMinObsInNode,
		    SEXP rdShrinkage,
		    SEXP rdBagFraction,
		    SEXP rcTrain,
		    SEXP rcTrainPatients,
		    SEXP rcFeatures):

		    dataConfig(radY,
		    radOffset,
		    radX,
		    raiXOrder,
		    rSorted,
		    rStrata,
		    radWeight,
		    radMisc,
		    rPriorCoeff,
		    rPatientId,
		    racVarClasses,
		    ralMonotoneVar,
		    rszFamily,
		    rdBagFraction,
		    rcTrain,
		    rcTrainPatients,
		    rcFeatures)
	{
		// Initialize DataDistParams
		InitszIRMeasure();

		// Initialize Tree Parameters
		treeConfig.cDepth = Rcpp::as<int>(rcDepth);
		treeConfig.cMinObsInNode = Rcpp::as<int>(rcMinObsInNode);
		treeConfig.dShrinkage = Rcpp::as<double>(rdShrinkage);
		treeConfig.cTrain = Rcpp::as<int>(rcTrain);

		Rcpp::NumericMatrix adX(radX);
		treeConfig.numColData =  adX.ncol();
	};


	//---------------------
	// Public destructor
	//---------------------
	~configStructs(){};

	//---------------------
	// Public Methods
	//---------------------
	DataDistParams& GetDataConfig() { return dataConfig;};
	TreeParams& GetTreeConfig() { return treeConfig;};

private:
	//-------------------
	// Private Variables
	//-------------------
	DataDistParams dataConfig;
	TreeParams treeConfig;

	//-------------------
	// Private Methods
	//-------------------
	inline void InitszIRMeasure()
	{
		// Check family specified
		if(dataConfig.family.empty())
		{
			throw GBM::failure("configStructs - Can't specify IR metric as family not initialized.");
		}
		// Get szIRMeasure and throw appropriate errors
		if(0 == dataConfig.family.compare(0, 8, "pairwise"))
		{
			std::size_t offsetToMeasure = dataConfig.family.find("_");
			if(offsetToMeasure == std::string::npos)
			{
				throw GBM::failure("Unable to locate IR metric required for pairwise");
			}
			const char* szIRMeasure = dataConfig.family.c_str() + offsetToMeasure + 1;
			dataConfig.szIRMeasure = szIRMeasure;
			dataConfig.family = "pairwise";
		}
		else
		{
			dataConfig.szIRMeasure = "";
		}


	}

};
#endif // CONFIGSTRUCTS_H
