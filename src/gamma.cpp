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
CGamma::CGamma(SEXP radMisc, const CDataset& data): CDistribution(radMisc, data)
{
}

//----------------------------------------
// Function Members - Public
//----------------------------------------
std::auto_ptr<CDistribution> CGamma::Create(SEXP radMisc, const CDataset& data,
										const char* szIRMeasure, int& cGroups, int& cTrain)
{

	return std::auto_ptr<CDistribution>(new CGamma(radMisc, data));
}

CGamma::~CGamma()
{
}


void CGamma::ComputeWorkingResponse
(
    const double *adF, 
    double *adZ, 
    const bag& afInBag,
    unsigned long nTrain
)
{
  unsigned long i = 0;
  double dF = 0.0;
  
  if (!(pData->y_ptr() && adF && adZ && pData->weight_ptr())) {
    throw GBM::invalid_argument();
  }

  for(i=0; i<nTrain; i++)
    {
      dF = adF[i] + ((pData->offset_ptr(false)==NULL) ? 0.0 : pData->offset_ptr(false)[i]);
      adZ[i] = pData->y_ptr()[i]*std::exp(-dF)-1.0;
    }
  
}


void CGamma::InitF
(
    double &dInitF, 
    unsigned long cLength
)
{
    double dSum=0.0;
    double dTotalWeight = 0.0;
    double Min = -19.0;
    double Max = +19.0;
    unsigned long i=0;

    if(pData->offset_ptr(false)==NULL)
    {
        for(i=0; i<cLength; i++)
        {
            dSum += pData->weight_ptr()[i]*pData->y_ptr()[i];
            dTotalWeight += pData->weight_ptr()[i];
        }
    }
    else
    {
        for(i=0; i<cLength; i++)
        {
	  dSum += pData->weight_ptr()[i]*pData->y_ptr()[i]*std::exp(-pData->offset_ptr(false)[i]);
	  dTotalWeight += pData->weight_ptr()[i];
        }
    }

    if (dSum <= 0.0) { dInitF = Min; }
    else { dInitF = std::log(dSum/dTotalWeight); }
    
    if (dInitF < Min) { dInitF = Min; }
    if (dInitF > Max) { dInitF = Max; }
}


double CGamma::Deviance
(
    const double *adF,
    unsigned long cLength,
    bool isValidationSet
)
{
  unsigned long i=0;
  double dL = 0.0;
  double dW = 0.0;
  double dF = 0.0;
  
  if(isValidationSet)
  {
	pData->shift_to_validation();
  }

  for(i=0; i!=cLength; i++)
    {
      dF = adF[i] + ((pData->offset_ptr(false)==NULL) ? 0.0 : pData->offset_ptr(false)[i]);
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
 const double *adF,
 double *adZ,
 const std::vector<unsigned long>& aiNodeAssign,
 unsigned long nTrain,
 VEC_P_NODETERMINAL vecpTermNodes,
 unsigned long cTermNodes,
 unsigned long cMinObsInNode,
 const bag& afInBag,
 const double *adFadj
)
{
 
  double dF = 0.0;
  unsigned long iObs = 0;
  unsigned long iNode = 0;
  double MaxVal = 19.0;
  double MinVal = -19.0;

  vecdNum.resize(cTermNodes);
  vecdNum.assign(vecdNum.size(),0.0);
  vecdDen.resize(cTermNodes);
  vecdDen.assign(vecdDen.size(),0.0);
  
  vecdMax.resize(cTermNodes);
  vecdMax.assign(vecdMax.size(),-HUGE_VAL);
  vecdMin.resize(cTermNodes);
  vecdMin.assign(vecdMin.size(),HUGE_VAL);
  
  for(iObs=0; iObs<nTrain; iObs++)
    {
      if(afInBag[iObs])
	{
	  dF = adF[iObs] + ((pData->offset_ptr(false)==NULL) ? 0.0 : pData->offset_ptr(false)[iObs]);
	  vecdNum[aiNodeAssign[iObs]] += pData->weight_ptr()[iObs]*pData->y_ptr()[iObs]*std::exp(-dF);
	  vecdDen[aiNodeAssign[iObs]] += pData->weight_ptr()[iObs];
	  
	  // Keep track of largest and smallest prediction in each node
	  vecdMax[aiNodeAssign[iObs]] = R::fmax2(dF,vecdMax[aiNodeAssign[iObs]]);
	  vecdMin[aiNodeAssign[iObs]] = R::fmin2(dF,vecdMin[aiNodeAssign[iObs]]);
	}
    }
  
  for(iNode=0; iNode<cTermNodes; iNode++)
    {
      if(vecpTermNodes[iNode]!=NULL)
	{
	  if(vecdNum[iNode] == 0.0)
	    {
	      // Taken from poisson.cpp
	      
	      // DEBUG: if vecdNum==0 then prediction = -Inf
	      // Not sure what else to do except plug in an arbitrary
	      //   negative number, -1? -10? Let's use -19, then make
	      //   sure |adF| < 19 always.
	      vecpTermNodes[iNode]->dPrediction = MinVal;
	    }
	  
	  else if(vecdDen[iNode] == 0.0) { vecpTermNodes[iNode]->dPrediction = 0.0; }
	  
	  else { vecpTermNodes[iNode]->dPrediction = std::log(vecdNum[iNode]/vecdDen[iNode]); }
	  
	  if (vecdMax[iNode]+vecpTermNodes[iNode]->dPrediction > MaxVal)
	    { vecpTermNodes[iNode]->dPrediction = MaxVal - vecdMax[iNode]; }
	  if (vecdMin[iNode]+vecpTermNodes[iNode]->dPrediction < MinVal)
	    { vecpTermNodes[iNode]->dPrediction = MinVal - vecdMin[iNode]; }
	  
	}
    }
}

double CGamma::BagImprovement
(
	const double *adF,
	const double *adFadj,
	const bag& afInBag,
	double dStepSize,
	unsigned long nTrain
)
{
	double dReturnValue = 0.0;
	double dF = 0.0;
	double dW = 0.0;
	unsigned long i = 0;

	for(i=0; i<nTrain; i++)
	{
		if(!afInBag[i])
		{
			dF = adF[i] + ((pData->offset_ptr(false)==NULL) ? 0.0 : pData->offset_ptr(false)[i]);
			dReturnValue += pData->weight_ptr()[i]*(pData->y_ptr()[i]*std::exp(-dF)*(1.0-exp(-dStepSize*adFadj[i])) - dStepSize*adFadj[i]);
			dW += pData->weight_ptr()[i];
		}
	}

	return 2*dReturnValue/dW;
}
