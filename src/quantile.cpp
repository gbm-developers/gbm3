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
CQuantile::CQuantile(double alpha): mpLocM("Other")
{
	dAlpha = alpha;
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
		throw GBM::failure("Quantile dist requires misc to initialization.");
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
	for(i=0; i<data.get_trainSize(); i++)
	{
		adZ[i] = (data.y_ptr()[i] > adF[i]+data.offset_ptr()[i]) ? dAlpha : -(1.0-dAlpha);
	}

}


double CQuantile::InitF
(
	const CDataset& data
)
{
    double dOffset=0.0;
    vecd.resize(data.get_trainSize());
    for(long i=0; i< data.get_trainSize(); i++)
    {
        dOffset = data.offset_ptr()[i];
        vecd[i] = data.y_ptr()[i] - dOffset;
    }

    return mpLocM.weightedQuantile(data.get_trainSize(), &vecd[0], data.weight_ptr(), dAlpha);
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
    long cLength = data.get_trainSize();
    if(isValidationSet)
    {
 	   data.shift_to_validation();
 	   cLength = data.GetValidSize();
    }


	for(i=0; i<cLength; i++)
	{
		if(data.y_ptr()[i] > adF[i] + data.offset_ptr()[i])
		{
			dL += data.weight_ptr()[i]*dAlpha*(data.y_ptr()[i] - adF[i]-data.offset_ptr()[i]);
		}
		else
		{
			dL += data.weight_ptr()[i]*(1.0-dAlpha)*(adF[i]+data.offset_ptr()[i] - data.y_ptr()[i]);
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

  vecd.resize(data.get_trainSize()); // should already be this size from InitF
  std::vector<double> adW2(data.get_trainSize());

  for(iNode=0; iNode<cTermNodes; iNode++)
    {
      if(treeComps.GetTermNodes()[iNode]->cN >= treeComps.GetMinNodeObs())
        {
	  iVecd = 0;
	  for(iObs=0; iObs< data.get_trainSize(); iObs++)
            {
	      if(data.GetBagElem(iObs) && (treeComps.GetNodeAssign()[iObs] == iNode))
                {
		  dOffset = data.offset_ptr()[iObs];
		  
		  vecd[iVecd] = data.y_ptr()[iObs] - dOffset - adF[iObs];
		  adW2[iVecd] = data.weight_ptr()[iObs];
		  iVecd++;
                }
            }
	  
	 treeComps.GetTermNodes()[iNode]->dPrediction = mpLocM.weightedQuantile(iVecd, &vecd[0], &adW2[0], dAlpha);
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

    for(i=0; i <data.get_trainSize(); i++)
    {
        if(!data.GetBagElem(i))
        {
            dF = adF[i] + data.offset_ptr()[i];

            if(data.y_ptr()[i] > dF)
            {
                dReturnValue += data.weight_ptr()[i]*dAlpha*(data.y_ptr()[i]-dF);
            }
            else
            {
                dReturnValue += data.weight_ptr()[i]*(1-dAlpha)*(dF-data.y_ptr()[i]);
            }

            if(data.y_ptr()[i] > dF+shrinkage*adFadj[i])
            {
                dReturnValue -= data.weight_ptr()[i]*dAlpha*
                                (data.y_ptr()[i] - dF-shrinkage*adFadj[i]);
            }
            else
            {
                dReturnValue -= data.weight_ptr()[i]*(1-dAlpha)*
                                (dF+shrinkage*adFadj[i] - data.y_ptr()[i]);
            }
            dW += data.weight_ptr()[i];
        }
    }
    return dReturnValue/dW;
}

