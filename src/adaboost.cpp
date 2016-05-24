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
CDistribution* CAdaBoost::Create(DataDistParams& distparams)
{
 	return new CAdaBoost();
}

CAdaBoost::~CAdaBoost()
{
}

void CAdaBoost::ComputeWorkingResponse
(
 const CDataset& data,
 const double* kFuncEstimate,
 double* residuals
)
{


	for(unsigned long i=0; i<data.get_trainsize(); i++)
	{
		residuals[i] = -(2*data.y_ptr()[i]-1) * std::exp(-(2*data.y_ptr()[i]-1)*(data.offset_ptr()[i]+kFuncEstimate[i]));
	}


}



double CAdaBoost::InitF
(
 const CDataset& data
)
{
    double numerator = 0.0;
    double denominator = 0.0;


	for(unsigned long i=0; i< data.get_trainsize(); i++)
	{
		if(data.y_ptr()[i]==1.0)
		{
			numerator += data.weight_ptr()[i] * std::exp(-data.offset_ptr()[i]);
		}
		else
		{
			denominator += data.weight_ptr()[i] * std::exp(data.offset_ptr()[i]);
		}
	}

    
    return 0.5*std::log(numerator/denominator);
}


double CAdaBoost::Deviance
(
	const CDataset& data,
    const double *adF,
    bool is_validationset
)
{
    unsigned long i=0;
    double loss = 0.0;
    double weight = 0.0;

    // Switch to validation set if necessary
    unsigned long num_of_rows_in_set = data.get_trainsize();
    if(is_validationset)
    {
    	data.shift_to_validation();
    	num_of_rows_in_set = data.get_validsize();
    }



	for(i=0; i!=num_of_rows_in_set; i++)
	{
		loss += data.weight_ptr()[i] * std::exp(-(2*data.y_ptr()[i]-1)*(data.offset_ptr()[i]+adF[i]));
		weight += data.weight_ptr()[i];
	}


    // Switch back to trainig set if necessary
   if(is_validationset)
   {
	   data.shift_to_train();
   }

   //TODO: Check if weights are all zero for validation set
   if((weight == 0.0) && (loss == 0.0))
   {
	   return nan("");
   }
   else if(weight == 0.0)
   {
	   return HUGE_VAL;
   }

   return loss/weight;
}


void CAdaBoost::FitBestConstant
(
	const CDataset& data,
    const double* adF,
    unsigned long cTermNodes,
    double* adZ,
    CTreeComps& treeComps
)
{
  double dF = 0.0;
  unsigned long iObs = 0;
  unsigned long iNode = 0;
  numerator_bestconstant_.resize(cTermNodes);
  numerator_bestconstant_.assign(numerator_bestconstant_.size(),0.0);
  denominator_bestconstant_.resize(cTermNodes);
  denominator_bestconstant_.assign(denominator_bestconstant_.size(),0.0);
    

  for(iObs=0; iObs< data.get_trainsize(); iObs++)
    {
      if(data.get_bag_element(iObs))
        {
	  dF = adF[iObs] + data.offset_ptr()[iObs];
	  numerator_bestconstant_[treeComps.get_node_assignments()[iObs]] +=
	    data.weight_ptr()[iObs]*(2*data.y_ptr()[iObs]-1)*std::exp(-(2*data.y_ptr()[iObs]-1)*dF);
	  denominator_bestconstant_[treeComps.get_node_assignments()[iObs]] +=
	    data.weight_ptr()[iObs]*std::exp(-(2*data.y_ptr()[iObs]-1)*dF);
        }
    }
  
  for(iNode=0; iNode<cTermNodes; iNode++)
    {
      if(treeComps.get_terminal_nodes()[iNode]!=NULL)
        {
	  if(denominator_bestconstant_[iNode] == 0)
            {
	      	  treeComps.get_terminal_nodes()[iNode]->prediction = 0.0;
            }
	  else
            {
	      treeComps.get_terminal_nodes()[iNode]->prediction =
		numerator_bestconstant_[iNode]/denominator_bestconstant_[iNode];
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
