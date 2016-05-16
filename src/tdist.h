//------------------------------------------------------------------------------
//
//  File:       tdist.h
//
//  Description:   Distribution object to implement t-distribution
//
//  History:    04/04/2008   Created
//
//------------------------------------------------------------------------------

#ifndef __tdistGBM_h__
#define __tdistGBM_h__

//------------------------------
// Includes
//------------------------------
#include "distribution.h"
#include "locationm.h"
#include <algorithm>
#include <memory>

//------------------------------
// Class definition
//------------------------------
class CTDist : public CDistribution
{
  
public:
	//---------------------
	// Factory Function
	//---------------------
	static CDistribution* Create(DataDistParams& distParams);

	//---------------------
	// Public destructor
	//---------------------
	virtual ~CTDist();
  
    //---------------------
    // Public Functions
    //---------------------
	void ComputeWorkingResponse(const CDataset& data,
			const double *adF,
			      	  double *adZ);

    double InitF(const CDataset& data);
    
    void FitBestConstant(const CDataset& data,
    		const double *adF,
			 unsigned long cTermNodes,
			 double* adZ,
			 CTreeComps& treeComps);

    double Deviance(const CDataset& data,
    				const double *adF,
                    bool isValidationSet=false);

    double BagImprovement(const CDataset& data,
    					  const double *adF,
    					  const bag& afInBag,
                          const double shrinkage,
                          const double* adFadj);

private:
    //----------------------
    // Private Constructors
    //----------------------
    CTDist(double nu);

	//-------------------
	// Private Variables
	//-------------------
    double mdNu;
    CLocationM mpLocM;
};

#endif // __tdistGBM_h__



