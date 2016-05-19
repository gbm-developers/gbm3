//------------------------------------------------------------------------------
//
//  File:       distribution.h
//
//  Description: Header file for distribution object used in GBM. This class
//    is the parent class from which all other distributions inherit.
//
//
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//				22/02/2016  jhickey: modified to implement factory pattern
//------------------------------------------------------------------------------

#ifndef __distribution_h__
#define __distribution_h__

//------------------------------
// Includes
//------------------------------
#include "configStructs.h"
#include "dataset.h"
#include "gbmTreeComps.h"
#include "node.h"
#include <vector>
#include <Rcpp.h>

//------------------------------
// Class definition
//------------------------------
class CDistribution
{

public:
	//----------------------
	// Public Constructors
	//----------------------
    CDistribution();

  	//---------------------
  	// Public destructor
  	//---------------------
    virtual ~CDistribution();

    //---------------------
    // Public Functions
    //---------------------
  	int GetNumGroups() const;
  	void SetNumGroups(int GroupVal);

  	// shifts the ptrs() as appropriate
  	template<typename T>
  	T* shift_ptr(T* x, std::ptrdiff_t y){
  		if(x)
  		{
  			return x+y;
  		}
  		else
  		{
  			return x;
  		}
  	}

     //---------------------
     // Public Virtual Functions
     //---------------------
     virtual void Initialize(const CDataset& data) { };
     virtual void ComputeWorkingResponse(const CDataset& data,
					const double *adF,
					double *adZ) = 0;

    virtual double InitF(const CDataset& data) = 0;

    virtual double Deviance(const CDataset& data, const double *adF,
                            bool isValidationSet=false) = 0;

    virtual void FitBestConstant(const CDataset& data, const double *adF,
						  unsigned long cTermNodes,
						  double* adZ, CTreeComps& treeComps) = 0;

    virtual double BagImprovement(const CDataset& data,
				  const double *adF,
				  const bag& afInBag,
				  const double shrinkage,
				  const double* adFadj) = 0;

    virtual void bagIt(CDataset& data);
 private:

    //---------------------
    // Private Variables
    //---------------------
    int cGroups;
};

#endif // __distribution_h__



