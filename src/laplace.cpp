//-----------------------------------
//
// File: laplace.cpp
//
// Description: laplace distribution for GBM.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "laplace.h"
#include <vector>

//----------------------------------------
// Function Members - Private
//----------------------------------------
CLaplace::CLaplace(): mpLocM("Other")
{
}

//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CLaplace::Create(DataDistParams& distParams)
{
	return new CLaplace();
}

CLaplace::~CLaplace()
{
}


void CLaplace::ComputeWorkingResponse
(
 const CDataset& data,
 const double *adF,
 double *adZ
)
{
    unsigned long i = 0;
	for(i=0; i<data.get_trainSize(); i++)
	{
		adZ[i] = (data.y_ptr()[i] - data.offset_ptr()[i] - adF[i]) > 0.0 ? 1.0 : -1.0;
	}

}



double CLaplace::InitF
(
	const CDataset& data
)
{
  double dOffset = 0.0;
  unsigned long ii = 0;

  std::vector<double> adArr(data.get_trainSize());
  
  for (ii = 0; ii < data.get_trainSize(); ii++)
    {
      dOffset = data.offset_ptr()[ii];
      adArr[ii] = data.y_ptr()[ii] - dOffset;
    }
  
  return mpLocM.weightedQuantile(data.get_trainSize(), &adArr[0], data.weight_ptr(), 0.5); // median
}


double CLaplace::Deviance
(
	const CDataset& data,
    const double *adF,
    bool isValidationSet
)
{
    unsigned long i=0;
    double dL = 0.0;
    double dW = 0.0;

    long cLength = data.get_trainSize();
    if(isValidationSet)
    {
    	data.shift_to_validation();
    	cLength = data.GetValidSize();
    }

    if(data.offset_ptr() == NULL)
    {
        for(i=0; i<cLength; i++)
        {
            dL += data.weight_ptr()[i]*fabs(data.y_ptr()[i]-adF[i]);
            dW += data.weight_ptr()[i];
        }
    }
    else
    {
        for(i=0; i<cLength; i++)
        {
            dL += data.weight_ptr()[i]*fabs(data.y_ptr()[i]-data.offset_ptr()[i]-adF[i]);
            dW += data.weight_ptr()[i];
        }
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


// DEBUG: needs weighted median
void CLaplace::FitBestConstant
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
  
//    vecd.resize(nTrain); // should already be this size from InitF

  std::vector<double> adArr(data.get_trainSize());
  std::vector<double> adW2(data.get_trainSize());
  
  for(iNode=0; iNode<cTermNodes; iNode++)
    {
      if(treeComps.GetTermNodes()[iNode]->cN >= treeComps.GetMinNodeObs())
        {
	  iVecd = 0;
	  for(iObs=0; iObs<data.get_trainSize(); iObs++)
            {
	      if(data.GetBagElem(iObs) && (treeComps.GetNodeAssign()[iObs] == iNode))
                {
		  dOffset =  data.offset_ptr()[iObs];
		  adArr[iVecd] = data.y_ptr()[iObs] - dOffset - adF[iObs];
		  adW2[iVecd] = data.weight_ptr()[iObs];
		  iVecd++;
		}
	      
            }
	  
	  treeComps.GetTermNodes()[iNode]->dPrediction = mpLocM.weightedQuantile(iVecd, &adArr[0], &adW2[0], 0.5); // median
	  
        }
    }
}



double CLaplace::BagImprovement
(
	const CDataset& data,
    const double *adF,
    const bag& afInBag,
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

            dReturnValue +=
                data.weight_ptr()[i]*(fabs(data.y_ptr()[i]-dF) - fabs(data.y_ptr()[i]-dF-shrinkage*adFadj[i]));
            dW += data.weight_ptr()[i];
        }
    }

    return dReturnValue/dW;
}
