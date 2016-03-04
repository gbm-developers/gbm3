//-----------------------------------
//
// File: gbmBuilder.cpp
//
// Description: class that builds the gbm_engine object
//
//-----------------------------------

//------------------------------
// Includes
//------------------------------
#include "gbmBuilder.h"

//----------------------------------------
// Function Members - Public
//----------------------------------------
//-----------------------------------
// Function: CGBMBuilder
//
// Returns: none
//
// Description: Default constructor for gbm builder - initializes pointer to GBM
//
// Parameters: none
//-----------------------------------
CGBMBuilder::CGBMBuilder()
{
	_gbmPtr = new CGBM();
	hasDataAndDist = false;
	hasTreeContainer = false;
}

//-----------------------------------
// Function: ~CGBMBuilder
//
// Returns: none
//
// Description: Default destructor for GBM builder
//
// Parameters: none
//-----------------------------------
CGBMBuilder::~CGBMBuilder()
{
}

//-----------------------------------
// Function: Build
//
// Returns: auto_ptr to built builder object
//
// Description: builds the gbm engine and returns a ptr to it.
//
// Parameters: none
//-----------------------------------
CGBM* CGBMBuilder::Build()
{
	if(hasDataAndDist && hasTreeContainer)
	{
		return _gbmPtr;
	}
	else
	{
		delete _gbmPtr;
		_gbmPtr = NULL;
		throw GBM::failure("GBM object could not be built - missing: "
				"data or distribution or weak learner");
	}
}

//-----------------------------------
// Function: BuildDataAndDist
//
// Returns: none
//
// Description: builds the data and distribution component of the gbm engine.
//
// Parameters:
//-----------------------------------
void CGBMBuilder::BuildDataAndDistribution(CDistribution* pDist)
{
	_gbmPtr-> SetDataAndDistribution(pDist);
	hasDataAndDist = true;
}

//-----------------------------------
// Function: BuildTreeContainer
//
// Returns: none
//
// Description: builds the tree components of the gbm engine.
//
// Parameters:
//-----------------------------------
void CGBMBuilder::BuildTreeContainer(double dLambda,
   	    unsigned long cTrain,
   	    unsigned long cFeatures,
   	    double dBagFraction,
   	    unsigned long cDepth,
   	    unsigned long cMinObsInNode,
   	    int cGroups)
{
	_gbmPtr-> SetTreeContainer(dLambda, cTrain, cFeatures,
    		dBagFraction, cDepth, cMinObsInNode, cGroups);
	hasTreeContainer = true;
}


