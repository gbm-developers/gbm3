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
#include "node.h"
#include "dataset.h"
#include "gbmTreeComps.h"
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
    CDistribution(SEXP radMisc);

  	//---------------------
  	// Public destructor
  	//---------------------
    virtual ~CDistribution();

    //---------------------
    // Public Functions
    //---------------------
  	bool has_misc() const ;
  	const double* misc_ptr(bool require=false) const;
  	double* misc_ptr(bool require=false);
  	int GetNumGroups() const;
  	void SetNumGroups(int GroupVal);

  	// shifts the misc_ptr() as appropriate
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
    virtual void Initialize(const CDataset* pData) { };

    virtual void ComputeWorkingResponse(const CDataset* pData,
    								const double *adF,
									double *adZ) = 0;

    virtual double InitF(const CDataset* pData) = 0;

    virtual double Deviance(const CDataset* pData, const double *adF,
                            bool isValidationSet=false) = 0;

    virtual void FitBestConstant(const CDataset* pData, const double *adF,
						  unsigned long cTermNodes,
						  double* adZ, CTreeComps* pTreeComps) = 0;

    virtual double BagImprovement(const CDataset& data,
    							  const double *adF,
								  const bag& afInBag,
								  const double shrinkage, const double* adFadj) = 0;

private:
    //---------------------
    // Private Variables
    //---------------------
    Rcpp::NumericVector adMisc;
    int cGroups;
    bool distHasMisc;
    bool distRequiresMisc;

};

#endif // __distribution_h__



