//-----------------------------------
//
// File: quantile.cpp
//
// Description: quantile distribution for GBM.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "quantile.h"

//----------------------------------------
// Function Members - Private
//----------------------------------------
CQuantile::CQuantile(double alpha): mplocm_("Other")
{
	alpha_ = alpha;
}


//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CQuantile::Create(DataDistParams& distParams)
{
	// Check that misc exists
	double alpha = Rcpp::as<double>(distParams.misc[0]);
	if(!GBM_FUNC::has_value(alpha))
	{
		throw GBM::Failure("Quantile dist requires misc to initialization.");
	}
	return new CQuantile(alpha);
}

CQuantile::~CQuantile()
{
}

void CQuantile::ComputeWorkingResponse
(
	const CDataset& data,
    const double *adF,
    double *adZ
)
{
    unsigned long i = 0;
	for(i=0; i<data.get_trainsize(); i++)
	{
		adZ[i] = (data.y_ptr()[i] > adF[i]+data.offset_ptr()[i]) ? alpha_ : -(1.0-alpha_);
	}

}


double CQuantile::InitF
(
	const CDataset& data
)
{
    double dOffset=0.0;
    vecd_.resize(data.get_trainsize());
    for(unsigned long i=0; i< data.get_trainsize(); i++)
    {
        dOffset = data.offset_ptr()[i];
        vecd_[i] = data.y_ptr()[i] - dOffset;
    }

    return mplocm_.weightedQuantile(data.get_trainsize(), &vecd_[0], data.weight_ptr(), alpha_);
}


double CQuantile::Deviance
(
	const CDataset& data,
    const double *adF,
    bool isValidationSet
)
{
    unsigned long i=0;
    double dL = 0.0;
    double dW = 0.0;

    // Switch to validation set if necessary
    unsigned long cLength = data.get_trainsize();
    if(isValidationSet)
    {
 	   data.shift_to_validation();
 	   cLength = data.get_validsize();
    }


	for(i=0; i<cLength; i++)
	{
		if(data.y_ptr()[i] > adF[i] + data.offset_ptr()[i])
		{
			dL += data.weight_ptr()[i]*alpha_*(data.y_ptr()[i] - adF[i]-data.offset_ptr()[i]);
		}
		else
		{
			dL += data.weight_ptr()[i]*(1.0-alpha_)*(adF[i]+data.offset_ptr()[i] - data.y_ptr()[i]);
		}
		dW += data.weight_ptr()[i];
	}


    // Switch back to training set if necessary
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

void CQuantile::FitBestConstant
(
	const CDataset& data,
    const double *adF,
    unsigned long cTermNodes,
    double* adZ,
    CTreeComps& treeComps
)
{
  unsigned long iNode = 0;
  unsigned long iObs = 0;
  unsigned long iVecd = 0;
  double dOffset;

  vecd_.resize(data.get_trainsize()); // should already be this size from InitF
  std::vector<double> adW2(data.get_trainsize());

  for(iNode=0; iNode<cTermNodes; iNode++)
    {
      if(treeComps.get_terminal_nodes()[iNode]->numobs >= treeComps.min_num_obs_required())
        {
	  iVecd = 0;
	  for(iObs=0; iObs< data.get_trainsize(); iObs++)
            {
	      if(data.get_bag_element(iObs) && (treeComps.get_node_assignments()[iObs] == iNode))
                {
		  dOffset = data.offset_ptr()[iObs];
		  
		  vecd_[iVecd] = data.y_ptr()[iObs] - dOffset - adF[iObs];
		  adW2[iVecd] = data.weight_ptr()[iObs];
		  iVecd++;
                }
            }
	  
	 treeComps.get_terminal_nodes()[iNode]->prediction = mplocm_.weightedQuantile(iVecd, &vecd_[0], &adW2[0], alpha_);
	}
    }
}


double CQuantile::BagImprovement
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

    for(i=0; i <data.get_trainsize(); i++)
    {
        if(!data.get_bag_element(i))
        {
            dF = adF[i] + data.offset_ptr()[i];

            if(data.y_ptr()[i] > dF)
            {
                dReturnValue += data.weight_ptr()[i]*alpha_*(data.y_ptr()[i]-dF);
            }
            else
            {
                dReturnValue += data.weight_ptr()[i]*(1-alpha_)*(dF-data.y_ptr()[i]);
            }

            if(data.y_ptr()[i] > dF+shrinkage*adFadj[i])
            {
                dReturnValue -= data.weight_ptr()[i]*alpha_*
                                (data.y_ptr()[i] - dF-shrinkage*adFadj[i]);
            }
            else
            {
                dReturnValue -= data.weight_ptr()[i]*(1-alpha_)*
                                (dF+shrinkage*adFadj[i] - data.y_ptr()[i]);
            }
            dW += data.weight_ptr()[i];
        }
    }
    return dReturnValue/dW;
}

