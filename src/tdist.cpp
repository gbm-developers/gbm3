//-----------------------------------
//
// File: tdist.cpp
//
// Description: t distribution implementation for GBM.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "locationm.h"
#include "tdist.h"
#include <vector>

//----------------------------------------
// Function Members - Private
//----------------------------------------
CTDist::CTDist(double nu):mpLocM("tdist", nu)
{
	mdNu = nu;
}


//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CTDist::Create(DataDistParams& distParams)
{
	// Check that misc exists
	double nu = Rcpp::as<double>(distParams.misc[0]);
	if(!GBM_FUNC::has_value(nu))
	{
		throw GBM::failure("T Dist requires misc to initialization.");
	}
	return new CTDist(nu);
}

CTDist::~CTDist()
{

}

void CTDist::ComputeWorkingResponse
(
	const CDataset& data,
    const double *adF,
    double *adZ
)
{
    unsigned long i = 0;
    double dU = 0.0;

	for(i=0; i<data.get_trainsize(); i++)
	{
		dU = data.y_ptr()[i] - data.offset_ptr()[i] - adF[i];
		adZ[i] = (2 * dU) / (mdNu + (dU * dU));
	}

}


double CTDist::InitF
(
	const CDataset& data
)
{

	// Get objects to pass into the LocM function
	std::vector<double> adArr(data.get_trainsize());

	for (unsigned long ii = 0; ii < data.get_trainsize(); ii++)
	{
		double dOffset = data.offset_ptr()[ii];
		adArr[ii] = data.y_ptr()[ii] - dOffset;
	}

	return mpLocM.LocationM(data.get_trainsize(), &adArr[0], data.weight_ptr(), 0.5);
}

double CTDist::Deviance
(
	const CDataset& data,
    const double *adF,
    bool isValidationSet
)
{
    unsigned long i=0;
    double dL = 0.0;
    double dW = 0.0;
	double dU = 0.0;

	// Switch to validation set if necessary
	unsigned long cLength = data.get_trainsize();
	if(isValidationSet)
	{
	   data.shift_to_validation();
	   cLength = data.get_validsize();
	}


	for(i=0; i<cLength; i++)
	{
		dU = data.y_ptr()[i] - data.offset_ptr()[i] - adF[i];
		dL += data.weight_ptr()[i] * std::log(mdNu + (dU * dU));
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


void CTDist::FitBestConstant
(
	const CDataset& data,
    const double *adF,
    unsigned long cTermNodes,
    double* adZ,
    CTreeComps& treeComps
)
{
   	// Local variables
    unsigned long iNode = 0;
    unsigned long iObs = 0;


    std::vector<double> adArr, adW;
	// Call LocM for the array of values on each node
    for(iNode=0; iNode<cTermNodes; iNode++)
    {
      if(treeComps.get_terminal_nodes()[iNode]->cN >= treeComps.min_num_obs_required())
        {
	  adArr.clear();
	  adW.clear();

	  for (iObs = 0; iObs < data.get_trainsize(); iObs++)
	    {
	      if(data.get_bag_element(iObs) && (treeComps.get_node_assignments()[iObs] == iNode))
                {
		  const double dOffset = data.offset_ptr()[iObs];
		  adArr.push_back(data.y_ptr()[iObs] - dOffset - adF[iObs]);
		  adW.push_back(data.weight_ptr()[iObs]);
                }
	    }

	  treeComps.get_terminal_nodes()[iNode]->dPrediction = mpLocM.LocationM(adArr.size(), &adArr[0],
							       &adW[0], 0.5);

        }
    }
}

double CTDist::BagImprovement
(
    const CDataset& data,
    const double *adF,
    const double shrinkage,
    const double* adFadj
)
{
    double dReturnValue = 0.0;
    unsigned long i = 0;
    double dW = 0.0;

    for(i=0; i<data.get_trainsize(); i++)
    {
        if(!data.get_bag_element(i))
        {
            const double dF = adF[i] + data.offset_ptr()[i];
	    const double dU = (data.y_ptr()[i] - dF);
	    const double dV = (data.y_ptr()[i] - dF - shrinkage * adFadj[i]) ;

            dReturnValue += data.weight_ptr()[i] * (std::log(mdNu + (dU * dU)) - log(mdNu + (dV * dV)));
            dW += data.weight_ptr()[i];
        }
    }

    return dReturnValue/dW;
}




