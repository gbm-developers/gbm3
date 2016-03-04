//------------------------------------------------------------------------------
//
//  File:       gbmBuilder.h
//
//  Description:   Header file for a class that builds the gbm object.
//
//------------------------------------------------------------------------------


#ifndef __gbmBuilder_h__
#define __gbmBuilder_h__

//------------------------------
// Includes
//------------------------------
#include "gbm_engine.h"
#include <vector>
#include <memory>

//------------------------------
// Class definition
//------------------------------
class CGBMBuilder
{
public:
	//----------------------
	// Public Constructors
	//----------------------
    CGBMBuilder();

	//---------------------
	// Public destructor
	//---------------------
    ~CGBMBuilder();

    //---------------------
	// Public Functions
	//---------------------
    CGBM* Build();
    void BuildDataAndDistribution(CDistribution* pDist);
    void BuildTreeContainer(double dLambda,
       	    unsigned long cTrain,
       	    unsigned long cFeatures,
       	    double dBagFraction,
       	    unsigned long cDepth,
       	    unsigned long cMinObsInNode,
       	    int cGroups);

private:
	//-------------------
	// Private Variables
	//-------------------
    CGBM* _gbmPtr;
    bool hasDataAndDist;
    bool hasTreeContainer;
};

#endif //  __gbmBuilder_h__
