//-----------------------------------
//
// File: adaboost.cpp
//
// Description: distribution used for adaboosting.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "adaboost.h"
#include <memory>

//----------------------------------------
// Function Members - Private
//----------------------------------------
CAdaBoost::CAdaBoost()
{
}

//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CAdaBoost::Create(DataDistParams& distParams)
{
 	return new CAdaBoost();
}

CAdaBoost::~CAdaBoost()
{
}

void CAdaBoost::ComputeWorkingResponse
(
 const CDataset& data,
 const double *adF,
 double *adZ
)
{


	for(unsigned long i=0; i<data.get_trainsize(); i++)
	{
		adZ[i] = -(2*data.y_ptr()[i]-1) * std::exp(-(2*data.y_ptr()[i]-1)*(data.offset_ptr()[i]+adF[i]));
	}


}



double CAdaBoost::InitF
(
 const CDataset& data
)
{
    double dNum = 0.0;
    double dDen = 0.0;


	for(unsigned long i=0; i< data.get_trainsize(); i++)
	{
		if(data.y_ptr()[i]==1.0)
		{
			dNum += data.weight_ptr()[i] * std::exp(-data.offset_ptr()[i]);
		}
		else
		{
			dDen += data.weight_ptr()[i] * std::exp(data.offset_ptr()[i]);
		}
	}

    
    return 0.5*std::log(dNum/dDen);
}


double CAdaBoost::Deviance
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



	for(i=0; i!=cLength; i++)
	{
		dL += data.weight_ptr()[i] * std::exp(-(2*data.y_ptr()[i]-1)*(data.offset_ptr()[i]+adF[i]));
		dW += data.weight_ptr()[i];
	}


    // Switch back to trainig set if necessary
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
	   return HUGE_VAL;
   }

   return dL/dW;
}


void CAdaBoost::FitBestConstant
(
	const CDataset& data,
    const double *adF,
    unsigned long cTermNodes,
    double* adZ,
    CTreeComps& treeComps
)
{
  double dF = 0.0;
  unsigned long iObs = 0;
  unsigned long iNode = 0;
  vecdNum.resize(cTermNodes);
  vecdNum.assign(vecdNum.size(),0.0);
  vecdDen.resize(cTermNodes);
  vecdDen.assign(vecdDen.size(),0.0);
    

  for(iObs=0; iObs< data.get_trainsize(); iObs++)
    {
      if(data.get_bag_element(iObs))
        {
	  dF = adF[iObs] + data.offset_ptr()[iObs];
	  vecdNum[treeComps.get_node_assignments()[iObs]] +=
	    data.weight_ptr()[iObs]*(2*data.y_ptr()[iObs]-1)*std::exp(-(2*data.y_ptr()[iObs]-1)*dF);
	  vecdDen[treeComps.get_node_assignments()[iObs]] +=
	    data.weight_ptr()[iObs]*std::exp(-(2*data.y_ptr()[iObs]-1)*dF);
        }
    }
  
  for(iNode=0; iNode<cTermNodes; iNode++)
    {
      if(treeComps.get_terminal_nodes()[iNode]!=NULL)
        {
	  if(vecdDen[iNode] == 0)
            {
	      	  treeComps.get_terminal_nodes()[iNode]->dPrediction = 0.0;
            }
	  else
            {
	      treeComps.get_terminal_nodes()[iNode]->dPrediction =
		vecdNum[iNode]/vecdDen[iNode];
            }
        }
    }
}


double CAdaBoost::BagImprovement
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

    for(i=0; i<data.get_trainsize(); i++)
    {
        if(!data.get_bag_element(i))
        {
            dF = adF[i] + data.offset_ptr()[i];

            dReturnValue += data.weight_ptr()[i]*
                (std::exp(-(2*data.y_ptr()[i]-1)*dF) -
                 std::exp(-(2*data.y_ptr()[i]-1)*(dF+(shrinkage)*(adFadj[i]))));
            dW += data.weight_ptr()[i];
        }
    }

    return dReturnValue/dW;
}
