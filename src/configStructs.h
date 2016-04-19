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

#ifndef __configStructs_h__
#define __configStructs_h__

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
		    SEXP radWeight,
		    SEXP radMisc,
		    SEXP racVarClasses,
		    SEXP ralMonotoneVar,
		    SEXP rszFamily,
		    SEXP rdBagFraction,
		    SEXP rcTrain,
		    SEXP rcFeatures): respY(radY)
{
		offset = radOffset;
		xValues = radX;
		xOrder = raiXOrder;
		varWeight =	radWeight;
		varClasses = racVarClasses;
		misc = radMisc;
		monotoneVar = ralMonotoneVar;
		cTrain = Rcpp::as<int>(rcTrain);
		cFeatures = Rcpp::as<int>(rcFeatures);
		dBagFraction = Rcpp::as<double>(rdBagFraction);
		family = Rcpp::as<std::string>(rszFamily);
		szIRMeasure = NULL;
}
	Rcpp::NumericMatrix respY;
	SEXP offset;
	SEXP xValues;
	SEXP xOrder;
	SEXP varWeight;
	SEXP varClasses;
	SEXP monotoneVar;
	SEXP misc;
	int cTrain;
	int cFeatures;
	double dBagFraction;
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
		    SEXP radWeight,
		    SEXP radMisc,
		    SEXP racVarClasses,
		    SEXP ralMonotoneVar,
		    SEXP rszFamily,
		    SEXP rcTrees,
		    SEXP rcDepth,
		    SEXP rcMinObsInNode,
		    SEXP rdShrinkage,
		    SEXP rdBagFraction,
		    SEXP rcTrain,
		    SEXP rcFeatures):
		    dataConfig(radY,
		    radOffset,
		    radX,
		    raiXOrder,
		    radWeight,
		    radMisc,
		    racVarClasses,
		    ralMonotoneVar,
		    rszFamily,
		    rdBagFraction,
		    rcTrain,
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
	DataDistParams GetDataConfig() const { return dataConfig;};
	TreeParams GetTreeConfig() const { return treeConfig;};

private:
	//-------------------
	// Private Variables
	//-------------------
	DataDistParams dataConfig;
	TreeParams treeConfig;

	//-------------------
	// Private Methods
	//-------------------
	void InitszIRMeasure()
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
#endif // __configStructs_h__*/
