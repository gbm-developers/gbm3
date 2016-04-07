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
#include <Rcpp.h>

//------------------------------
// Struct definitions
//------------------------------
struct DataDistParams
{
	SEXP respY;
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
		    SEXP rcFeatures)
	{
		// Initialize DataDistParams
		dataConfig.respY = radY;
		dataConfig.offset = radOffset;
		dataConfig.xValues = radX;
		dataConfig.xOrder = raiXOrder;
		dataConfig.varWeight =	radWeight;
		dataConfig.varClasses = racVarClasses;
		dataConfig.misc = radMisc;
		dataConfig.monotoneVar = ralMonotoneVar;
		dataConfig.cTrain = Rcpp::as<int>(rcTrain);
		dataConfig.cFeatures = Rcpp::as<int>(rcFeatures);
		dataConfig.dBagFraction = Rcpp::as<double>(rdBagFraction);
		dataConfig.family = Rcpp::as<std::string>(rszFamily);


		// Initialize Tree Parameter
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

};
#endif // __configStructs_h__*/
