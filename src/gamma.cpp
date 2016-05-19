//-----------------------------------
//
// File: gamma.cpp
//
// Description: Gamma distribution with natural
//  log link function ( mean = std::exp(prediction) )
//
// Notes: -19 <= prediction <= +19
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------

#include "gamma.h"
#include <math.h>
#include <iostream>

//----------------------------------------
// Function Members - Private
//----------------------------------------
CGamma::CGamma()
{
}

//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CGamma::Create(DataDistParams& distParams)
{

	return new CGamma();
}

CGamma::~CGamma()
{
}


void CGamma::ComputeWorkingResponse
(
	const CDataset& data,
    const double *adF, 
    double *adZ
)
{
  unsigned long i = 0;
  double dF = 0.0;
  
  if (!(data.y_ptr() && adF && adZ && data.weight_ptr())) {
    throw GBM::invalid_argument();
  }

  for(i=0; i<data.get_trainSize(); i++)
    {
      dF = adF[i] + data.offset_ptr()[i];
      adZ[i] = data.y_ptr()[i]*std::exp(-dF)-1.0;
    }
  
}


double CGamma::InitF
(
	const CDataset& data
)
{
    double dSum=0.0;
    double dTotalWeight = 0.0;
    double Min = -19.0;
    double Max = +19.0;
    unsigned long i=0;
    double dInitF = 0.0;

	for(i=0; i<data.get_trainSize(); i++)
	{
		dSum += data.weight_ptr()[i]*data.y_ptr()[i]*std::exp(-data.offset_ptr()[i]);
		dTotalWeight += data.weight_ptr()[i];
	}


    if (dSum <= 0.0) { dInitF = Min; }
    else { dInitF = std::log(dSum/dTotalWeight); }
    
    if (dInitF < Min) { dInitF = Min; }
    if (dInitF > Max) { dInitF = Max; }
    return dInitF;
}


double CGamma::Deviance
(
	const CDataset& data,
    const double *adF,
    bool isValidationSet
)
{
  unsigned long i=0;
  double dL = 0.0;
  double dW = 0.0;
  double dF = 0.0;
  
  long cLength = data.get_trainSize();
  if(isValidationSet)
  {
	data.shift_to_validation();
	cLength = data.GetValidSize();
  }

  for(i=0; i!=cLength; i++)
    {
      dF = adF[i] + data.offset_ptr()[i];
      dL += data.weight_ptr()[i]*(data.y_ptr()[i]*std::exp(-dF) + dF);
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

  return 2*dL/dW;
}


void CGamma::FitBestConstant
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
  double MaxVal = 19.0;
  double MinVal = -19.0;

  vector<double> vecdNum(cTermNodes, 0.0);
  vector<double> vecdDen(cTermNodes, 0.0);
  vector<double> vecdMax(cTermNodes, -HUGE_VAL);
  vector<double> vecdMin(cTermNodes, HUGE_VAL);

  for(iObs=0; iObs<data.get_trainSize(); iObs++)
    {
      if(data.GetBagElem(iObs))
	{
	  dF = adF[iObs] + data.offset_ptr()[iObs];
	  vecdNum[treeComps.GetNodeAssign()[iObs]] += data.weight_ptr()[iObs]*data.y_ptr()[iObs]*std::exp(-dF);
	  vecdDen[treeComps.GetNodeAssign()[iObs]] += data.weight_ptr()[iObs];
	  
	  // Keep track of largest and smallest prediction in each node
	  vecdMax[treeComps.GetNodeAssign()[iObs]] = R::fmax2(dF,vecdMax[treeComps.GetNodeAssign()[iObs]]);
	  vecdMin[treeComps.GetNodeAssign()[iObs]] = R::fmin2(dF,vecdMin[treeComps.GetNodeAssign()[iObs]]);
	}
    }
  
  for(iNode=0; iNode<cTermNodes; iNode++)
    {
      if(treeComps.GetTermNodes()[iNode]!=NULL)
	{
	  if(vecdNum[iNode] == 0.0)
	    {
	      // Taken from poisson.cpp
	      
	      // DEBUG: if vecdNum==0 then prediction = -Inf
	      // Not sure what else to do except plug in an arbitrary
	      //   negative number, -1? -10? Let's use -19, then make
	      //   sure |adF| < 19 always.
		  treeComps.GetTermNodes()[iNode]->dPrediction = MinVal;
	    }
	  
	  else if(vecdDen[iNode] == 0.0) { treeComps.GetTermNodes()[iNode]->dPrediction = 0.0; }
	  
	  else { treeComps.GetTermNodes()[iNode]->dPrediction = std::log(vecdNum[iNode]/vecdDen[iNode]); }
	  
	  if (vecdMax[iNode]+treeComps.GetTermNodes()[iNode]->dPrediction > MaxVal)
	    { treeComps.GetTermNodes()[iNode]->dPrediction = MaxVal - vecdMax[iNode]; }
	  if (vecdMin[iNode]+treeComps.GetTermNodes()[iNode]->dPrediction < MinVal)
	    { treeComps.GetTermNodes()[iNode]->dPrediction = MinVal - vecdMin[iNode]; }
	  
	}
    }
}

double CGamma::BagImprovement
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
			dF = adF[i] +  data.offset_ptr()[i];
			dReturnValue += data.weight_ptr()[i]*(data.y_ptr()[i]*std::exp(-dF)*(1.0-exp(-shrinkage*adFadj[i])) - shrinkage*adFadj[i]);
			dW += data.weight_ptr()[i];
		}
	}

	return 2*dReturnValue/dW;
}
