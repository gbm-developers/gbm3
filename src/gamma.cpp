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
CGamma::CGamma(SEXP radMisc): CDistribution(radMisc)
{
}

//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CGamma::Create(SEXP radMisc,
							  const char* szIRMeasure,  int& cTrain)
{

	return new CGamma(radMisc);
}

CGamma::~CGamma()
{
}


void CGamma::ComputeWorkingResponse
(
	const CDataset* pData,
    const double *adF, 
    double *adZ
)
{
  unsigned long i = 0;
  double dF = 0.0;
  
  if (!(pData->y_ptr() && adF && adZ && pData->weight_ptr())) {
    throw GBM::invalid_argument();
  }

  for(i=0; i<pData->get_trainSize(); i++)
    {
      dF = adF[i] + pData->offset_ptr(false)[i];
      adZ[i] = pData->y_ptr()[i]*std::exp(-dF)-1.0;
    }
  
}


double CGamma::InitF
(
	const CDataset* pData
)
{
    double dSum=0.0;
    double dTotalWeight = 0.0;
    double Min = -19.0;
    double Max = +19.0;
    unsigned long i=0;
    double dInitF = 0.0;

	for(i=0; i<pData->get_trainSize(); i++)
	{
		dSum += pData->weight_ptr()[i]*pData->y_ptr()[i]*std::exp(-pData->offset_ptr(false)[i]);
		dTotalWeight += pData->weight_ptr()[i];
	}


    if (dSum <= 0.0) { dInitF = Min; }
    else { dInitF = std::log(dSum/dTotalWeight); }
    
    if (dInitF < Min) { dInitF = Min; }
    if (dInitF > Max) { dInitF = Max; }
    return dInitF;
}


double CGamma::Deviance
(
	const CDataset* pData,
    const double *adF,
    bool isValidationSet
)
{
  unsigned long i=0;
  double dL = 0.0;
  double dW = 0.0;
  double dF = 0.0;
  
  long cLength = pData->get_trainSize();
  if(isValidationSet)
  {
	pData->shift_to_validation();
	cLength = pData->GetValidSize();
  }

  for(i=0; i!=cLength; i++)
    {
      dF = adF[i] + pData->offset_ptr(false)[i];
      dL += pData->weight_ptr()[i]*(pData->y_ptr()[i]*std::exp(-dF) + dF);
      dW += pData->weight_ptr()[i];
    }
  
  if(isValidationSet)
  {
	pData->shift_to_train();
  }

  return 2*dL/dW;
}


void CGamma::FitBestConstant
(
const CDataset* pData,
 const double *adF,
 unsigned long cTermNodes,
 double* adZ,
 CTreeComps* pTreeComps
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

  for(iObs=0; iObs<pData->get_trainSize(); iObs++)
    {
      if(pData->GetBagElem(iObs))
	{
	  dF = adF[iObs] + pData->offset_ptr(false)[iObs];
	  vecdNum[pTreeComps->GetNodeAssign()[iObs]] += pData->weight_ptr()[iObs]*pData->y_ptr()[iObs]*std::exp(-dF);
	  vecdDen[pTreeComps->GetNodeAssign()[iObs]] += pData->weight_ptr()[iObs];
	  
	  // Keep track of largest and smallest prediction in each node
	  vecdMax[pTreeComps->GetNodeAssign()[iObs]] = R::fmax2(dF,vecdMax[pTreeComps->GetNodeAssign()[iObs]]);
	  vecdMin[pTreeComps->GetNodeAssign()[iObs]] = R::fmin2(dF,vecdMin[pTreeComps->GetNodeAssign()[iObs]]);
	}
    }
  
  for(iNode=0; iNode<cTermNodes; iNode++)
    {
      if(pTreeComps->GetTermNodes()[iNode]!=NULL)
	{
	  if(vecdNum[iNode] == 0.0)
	    {
	      // Taken from poisson.cpp
	      
	      // DEBUG: if vecdNum==0 then prediction = -Inf
	      // Not sure what else to do except plug in an arbitrary
	      //   negative number, -1? -10? Let's use -19, then make
	      //   sure |adF| < 19 always.
		  pTreeComps->GetTermNodes()[iNode]->dPrediction = MinVal;
	    }
	  
	  else if(vecdDen[iNode] == 0.0) { pTreeComps->GetTermNodes()[iNode]->dPrediction = 0.0; }
	  
	  else { pTreeComps->GetTermNodes()[iNode]->dPrediction = std::log(vecdNum[iNode]/vecdDen[iNode]); }
	  
	  if (vecdMax[iNode]+pTreeComps->GetTermNodes()[iNode]->dPrediction > MaxVal)
	    { pTreeComps->GetTermNodes()[iNode]->dPrediction = MaxVal - vecdMax[iNode]; }
	  if (vecdMin[iNode]+pTreeComps->GetTermNodes()[iNode]->dPrediction < MinVal)
	    { pTreeComps->GetTermNodes()[iNode]->dPrediction = MinVal - vecdMin[iNode]; }
	  
	}
    }
}

double CGamma::BagImprovement
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
			dF = adF[i] +  data.offset_ptr(false)[i];
			dReturnValue += data.weight_ptr()[i]*(data.y_ptr()[i]*std::exp(-dF)*(1.0-exp(-shrinkage*adFadj[i])) - shrinkage*adFadj[i]);
			dW += data.weight_ptr()[i];
		}
	}

	return 2*dReturnValue/dW;
}
