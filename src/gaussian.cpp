//-----------------------------------
//
// File: gaussian.cpp
//
// Description: gaussian distribution implementation for GBM.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "gaussian.h"

//----------------------------------------
// Function Members - Private
//----------------------------------------
CGaussian::CGaussian()
{

}

//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CGaussian::Create(DataDistParams& distParams)
{
	return new CGaussian();
}

CGaussian::~CGaussian()
{
}


void CGaussian::ComputeWorkingResponse
(
 const CDataset& data,
 const double *adF,
 double *adZ
 )
{
  unsigned long i = 0;
  
  if (!(data.y_ptr() && adF && adZ && data.weight_ptr())) {
    throw GBM::invalid_argument();
  }
  

	for(i=0; i<data.get_trainSize(); i++)
	{
		adZ[i] = data.y_ptr()[i] - data.offset_ptr()[i] - adF[i];
	}

}

double CGaussian::InitF
(
	const CDataset& data
)
{
    double dSum=0.0;
    double dTotalWeight = 0.0;
    unsigned long i=0;

    // compute the mean

	for(i=0; i<data.get_trainSize(); i++)
	{
		dSum += data.weight_ptr()[i]*(data.y_ptr()[i] - data.offset_ptr()[i]);
		dTotalWeight += data.weight_ptr()[i];
	}


    return dSum/dTotalWeight;
}


double CGaussian::Deviance
(
	const CDataset& data,
    const double *adF,
    bool isValidationSet
)
{
    unsigned long i=0;
    double dL = 0.0;
    double dW = 0.0;

    unsigned long cLength = data.get_trainSize();
    if(isValidationSet)
    {
    	data.shift_to_validation();
    	cLength = data.GetValidSize();
    }



	for(i=0; i<cLength; i++)
	{
		dL += data.weight_ptr()[i]*(data.y_ptr()[i]-data.offset_ptr()[i]-adF[i])*
						  (data.y_ptr()[i]-data.offset_ptr()[i]-adF[i]);
		dW += data.weight_ptr()[i];
	}


    if(isValidationSet)
    {
    	data.shift_to_train();
    }

    //TODO: Check if weights are all zero for validation set
   if((dW == 0.0) && (dL == 0.0))
   {
	   return nan("");
   }
   else if(dW == 0.0)
   {
	   return copysign(HUGE_VAL, dL);
   }

    return dL/dW;
}


void CGaussian::FitBestConstant
(
	const CDataset& data,
    const double *adF,
    unsigned long cTermNodes,
    double* adZ,
    CTreeComps& treeComps
)
{
  // the tree aready stores the mean prediction
  // no refitting necessary
}

double CGaussian::BagImprovement
(
    const CDataset& data,
    const double *adF,
    const double shrinkage,
    const double* adFadj
)
{
    double dReturnValue = 0.0;
    double dF = 0.0;
    double dW = 0.0;
    unsigned long i = 0;

    for(i=0; i<data.get_trainSize(); i++)
    {
        if(!data.GetBagElem(i))
        {
            dF = adF[i] + data.offset_ptr()[i];

            dReturnValue += data.weight_ptr()[i]*shrinkage*adFadj[i]*
                            (2.0*(data.y_ptr()[i]-dF) - shrinkage*adFadj[i]);
            dW += data.weight_ptr()[i];
        }
    }

    return dReturnValue/dW;
}



