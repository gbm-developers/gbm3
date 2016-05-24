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
CDistribution* CPoisson::Create(DataDistParams& distParams)
{
	return new CPoisson();
}


CPoisson::~CPoisson()
{
}


void CPoisson::ComputeWorkingResponse
(
	const CDataset& data,
    const double *adF,
    double *adZ
)
{
    unsigned long i = 0;
    double dF = 0.0;

    // compute working response
    for(i=0; i < data.get_trainsize(); i++)
    {
        dF = adF[i] +  data.offset_ptr()[i];
        adZ[i] = data.y_ptr()[i] - std::exp(dF);
    }
}



double CPoisson::InitF
(
	const CDataset& data
)
{
    double dSum = 0.0;
    double dDenom = 0.0;
    unsigned long i = 0;


	for(i=0; i<data.get_trainsize(); i++)
	{
		dSum += data.weight_ptr()[i]*data.y_ptr()[i];
		dDenom += data.weight_ptr()[i]*std::exp(data.offset_ptr()[i]);
	}

    return std::log(dSum/dDenom);
}


double CPoisson::Deviance
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
		dL += data.weight_ptr()[i]*(data.y_ptr()[i]*(data.offset_ptr()[i]+adF[i]) -
						   std::exp(data.offset_ptr()[i]+adF[i]));
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
	   return copysign(HUGE_VAL, -dL);
   }
    return -2*dL/dW;
}


void CPoisson::FitBestConstant
(
	const CDataset& data,
    const double *adF,
    unsigned long cTermNodes,
    double* adZ,
    CTreeComps& treeComps
)
{
    unsigned long iObs = 0;
    unsigned long iNode = 0;
    vector<double> vecdNum(cTermNodes, 0.0);
    vector<double> vecdDen(cTermNodes, 0.0);
    vector<double> vecdMax(cTermNodes, -HUGE_VAL);
    vector<double> vecdMin(cTermNodes, HUGE_VAL);

	for(iObs=0; iObs<data.get_trainsize(); iObs++)
	{
		if(data.get_bag_element(iObs))
		{
			vecdNum[treeComps.get_node_assignments()[iObs]] += data.weight_ptr()[iObs]*data.y_ptr()[iObs];
			vecdDen[treeComps.get_node_assignments()[iObs]] +=
				data.weight_ptr()[iObs]*std::exp(data.offset_ptr()[iObs]+adF[iObs]);
		}
	}

    for(iNode=0; iNode<cTermNodes; iNode++)
    {
        if(treeComps.get_terminal_nodes()[iNode]!=NULL)
        {
            if(vecdNum[iNode] == 0.0)
            {
                // DEBUG: if vecdNum==0 then prediction = -Inf
                // Not sure what else to do except plug in an arbitrary
                //   negative number, -1? -10? Let's use -1, then make
                //   sure |adF| < 19 always.
            	treeComps.get_terminal_nodes()[iNode]->prediction = -19.0;
            }
            else if(vecdDen[iNode] == 0.0)
            {
            	treeComps.get_terminal_nodes()[iNode]->prediction = 0.0;
            }
            else
            {
            	treeComps.get_terminal_nodes()[iNode]->prediction =
                    std::log(vecdNum[iNode]/vecdDen[iNode]);
            }
            treeComps.get_terminal_nodes()[iNode]->prediction =
               R::fmin2(treeComps.get_terminal_nodes()[iNode]->prediction,
                     19-vecdMax[iNode]);
            treeComps.get_terminal_nodes()[iNode]->prediction =
               R::fmax2(treeComps.get_terminal_nodes()[iNode]->prediction,
                     -19-vecdMin[iNode]);
        }
    }
}


double CPoisson::BagImprovement
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
                            (data.y_ptr()[i]*shrinkage*adFadj[i] -
                             std::exp(dF+shrinkage*adFadj[i]) +
                             std::exp(dF));
            dW += data.weight_ptr()[i];
        }
    }

    return dReturnValue/dW;
}


