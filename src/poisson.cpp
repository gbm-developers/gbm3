//-----------------------------------
//
// File: poisson.cpp
//
// Description: poisson distribution for GBM.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "poisson.h"

//----------------------------------------
// Function Members - Private
//----------------------------------------
CPoisson::CPoisson()
{
}


//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CPoisson::Create(DataDistParams& distparams)
{
	return new CPoisson();
}


CPoisson::~CPoisson()
{
}


void CPoisson::ComputeWorkingResponse
(
	const CDataset& kData,
    const double* kFuncEstimate,
    double* residuals
)
{
    unsigned long i = 0;
    double delta_func_est = 0.0;

    // compute working response
    for(i=0; i < kData.get_trainsize(); i++)
    {
        delta_func_est = kFuncEstimate[i] +  kData.offset_ptr()[i];
        residuals[i] = kData.y_ptr()[i] - std::exp(delta_func_est);
    }
}



double CPoisson::InitF
(
	const CDataset& kData
)
{
    double sum = 0.0;
    double denom = 0.0;
    unsigned long i = 0;


	for(i=0; i<kData.get_trainsize(); i++)
	{
		sum += kData.weight_ptr()[i]*kData.y_ptr()[i];
		denom += kData.weight_ptr()[i]*std::exp(kData.offset_ptr()[i]);
	}

    return std::log(sum/denom);
}


double CPoisson::Deviance
(
	const CDataset& kData,
    const double* kFuncEstimate,
    bool is_validationset
)
{
    unsigned long i=0;
    double loss = 0.0;
    double weight = 0.0;

    // Switch to validation set if necessary
    unsigned long num_rows_in_set = kData.get_trainsize();
    if(is_validationset)
    {
 	   kData.shift_to_validation();
 	   num_rows_in_set = kData.get_validsize();
    }


	for(i=0; i<num_rows_in_set; i++)
	{
		loss += kData.weight_ptr()[i]*(kData.y_ptr()[i]*(kData.offset_ptr()[i]+kFuncEstimate[i]) -
						   std::exp(kData.offset_ptr()[i]+kFuncEstimate[i]));
		weight += kData.weight_ptr()[i];
   }


    // Switch back to training set if necessary
    if(is_validationset)
    {
 	   kData.shift_to_train();
    }

    //TODO: Check if weights are all zero for validation set
   if((weight == 0.0) && (loss == 0.0))
   {
	   return nan("");
   }
   else if(weight == 0.0)
   {
	   return copysign(HUGE_VAL, -loss);
   }
    return -2*loss/weight;
}


void CPoisson::FitBestConstant
(
	const CDataset& kData,
    const double* kFuncEstimate,
    unsigned long num_terminalnodes,
    double* residuals,
    CCARTTree& tree
)
{
    unsigned long obs_num = 0;
    unsigned long node_num = 0;
    vector<double> numerator_vec(num_terminalnodes, 0.0);
    vector<double> denominator_vec(num_terminalnodes, 0.0);
    vector<double> max_vec(num_terminalnodes, -HUGE_VAL);
    vector<double> min_vec(num_terminalnodes, HUGE_VAL);

	for(obs_num=0; obs_num<kData.get_trainsize(); obs_num++)
	{
		if(kData.get_bag_element(obs_num))
		{
			numerator_vec[tree.get_node_assignments()[obs_num]] += kData.weight_ptr()[obs_num]*kData.y_ptr()[obs_num];
			denominator_vec[tree.get_node_assignments()[obs_num]] +=
				kData.weight_ptr()[obs_num]*std::exp(kData.offset_ptr()[obs_num]+kFuncEstimate[obs_num]);
		}
	}

    for(node_num=0; node_num<num_terminalnodes; node_num++)
    {
        if(tree.get_terminal_nodes()[node_num]!=NULL)
        {
            if(numerator_vec[node_num] == 0.0)
            {
                // DEBUG: if vecdNum==0 then prediction = -Inf
                // Not sure what else to do except plug in an arbitrary
                //   negative number, -1? -10? Let's use -1, then make
                //   sure |adF| < 19 always.
            	tree.get_terminal_nodes()[node_num]->prediction = -19.0;
            }
            else if(denominator_vec[node_num] == 0.0)
            {
            	tree.get_terminal_nodes()[node_num]->prediction = 0.0;
            }
            else
            {
            	tree.get_terminal_nodes()[node_num]->prediction =
                    std::log(numerator_vec[node_num]/denominator_vec[node_num]);
            }
            tree.get_terminal_nodes()[node_num]->prediction =
               R::fmin2(tree.get_terminal_nodes()[node_num]->prediction,
                     19-max_vec[node_num]);
            tree.get_terminal_nodes()[node_num]->prediction =
               R::fmax2(tree.get_terminal_nodes()[node_num]->prediction,
                     -19-min_vec[node_num]);
        }
    }
}


double CPoisson::BagImprovement
(
	const CDataset& kData,
    const double* kFuncEstimate,
    const double kShrinkage,
    const double* kDeltaEstimate
)
{
    double returnvalue = 0.0;
    double delta_func_est = 0.0;
    double weight = 0.0;
    unsigned long i = 0;

    for(i=0; i<kData.get_trainsize(); i++)
    {
        if(!kData.get_bag_element(i))
        {
            delta_func_est = kFuncEstimate[i] + kData.offset_ptr()[i];

            returnvalue += kData.weight_ptr()[i]*
                            (kData.y_ptr()[i]*kShrinkage*kDeltaEstimate[i] -
                             std::exp(delta_func_est+kShrinkage*kDeltaEstimate[i]) +
                             std::exp(delta_func_est));
            weight += kData.weight_ptr()[i];
        }
    }

    return returnvalue/weight;
}


