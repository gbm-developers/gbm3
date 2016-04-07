//-----------------------------------
//
// File: huberized.cpp
//
// Description: huberized hinge loss.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "huberized.h"

//----------------------------------------
// Function Members - Private
//----------------------------------------
CHuberized::CHuberized(SEXP radMisc): CDistribution(radMisc)
{
}

//----------------------------------------
// Function Members - Public
//----------------------------------------

CDistribution* CHuberized::Create(SEXP radMisc,
								 const char* szIRMeasure,
								 int& cTrain)
{
	return new CHuberized(radMisc);
}

CHuberized::~CHuberized()
{
}


void CHuberized::ComputeWorkingResponse
(
	const CDataset* pData,
    const double *adF,
    double *adZ
)
{
   unsigned long i = 0;
   double dF = 0.0;

   for(i=0; i<pData->get_trainSize(); i++)
   {
      dF = adF[i] + pData->offset_ptr(false)[i];
      if( (2*pData->y_ptr()[i]-1)*dF < -1)
      {
         adZ[i] = -4 * (2*pData->y_ptr()[i]-1);
      }
      else if ( 1 - (2*pData->y_ptr()[i]-1)*dF < 0 )
      {
         adZ[i] = 0;
      }
      else
      {
         adZ[i] = -2 * (2*pData->y_ptr()[i]-1) * ( 1 - (2*pData->y_ptr()[i]-1)*dF );
      }
   }
}

double CHuberized::InitF
(
	const CDataset* pData
)
{
    unsigned long i=0;
    double dNum = 0.0;
    double dDen = 0.0;

    for(i=0; i<pData->get_trainSize(); i++)
    {
        if(pData->y_ptr()[i]==1.0)
        {
            dNum += pData->weight_ptr()[i];
        }
        else
        {
            dDen += pData->weight_ptr()[i];
        }
    }

    return dNum/dDen;
}


double CHuberized::Deviance
(
	const CDataset* pData,
    const double *adF,
    bool isValidationSet
)
{
   unsigned long i=0;
   double dL = 0.0;
   double dF = 0.0;
   double dW = 0.0;

   long cLength = pData->get_trainSize();
   if(isValidationSet)
   {
	   pData->shift_to_validation();
	   cLength = pData->GetValidSize();
   }


  for(i=0; i<cLength; i++)
  {
	 dF = pData->offset_ptr(false)[i]+adF[i];
	 if ( (2*pData->y_ptr()[i]-1)*adF[i] < -1 )
	 {
		dL += -pData->weight_ptr()[i]*4*(2*pData->y_ptr()[i]-1)*dF;
		dW += pData->weight_ptr()[i];
	 }
	 else if ( 1 - (2*pData->y_ptr()[i]-1)*dF < 0 )
	 {
		dL += 0;
		dW += pData->weight_ptr()[i];
	 }
	 else
	 {
		dL += pData->weight_ptr()[i] * ( 1 - (2*pData->y_ptr()[i]-1)*dF ) *
							( 1 - (2*pData->y_ptr()[i]-1)*dF );
		dW += pData->weight_ptr()[i];
	 }
  } // close for(


   if(isValidationSet)
   {
	   pData->shift_to_train();
   }

   return dL/dW;
}


void CHuberized::FitBestConstant
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
  
  vector<double> vecdNum(cTermNodes, 0.0);
  vector<double> vecdDen(cTermNodes, 0.0);

  for(iObs=0; iObs<pData->get_trainSize(); iObs++)
    {
      if(pData->GetBagElem(iObs))
        {
	  dF = adF[iObs] +  pData->offset_ptr(false)[iObs];
	  if( (2*pData->y_ptr()[iObs]-1)*adF[iObs] < -1 )
	  {
	    vecdNum[pTreeComps->GetNodeAssign()[iObs]] +=
	      pData->weight_ptr()[iObs]*4*(2*pData->y_ptr()[iObs]-1);
	    vecdDen[pTreeComps->GetNodeAssign()[iObs]] +=
	      -pData->weight_ptr()[iObs]*4*(2*pData->y_ptr()[iObs]-1)*dF;
	  }
	  else if ( 1 - (2*pData->y_ptr()[iObs]-1)*adF[iObs] < 0 ){
	    vecdNum[pTreeComps->GetNodeAssign()[iObs]] += 0;
	    vecdDen[pTreeComps->GetNodeAssign()[iObs]] += 0;
	  }
	  else{
	    vecdNum[pTreeComps->GetNodeAssign()[iObs]] += pData->weight_ptr()[iObs]*2*(2*pData->y_ptr()[iObs]-1)*( 1 - (2*pData->y_ptr()[iObs]-1)*adF[iObs] );
	    vecdDen[pTreeComps->GetNodeAssign()[iObs]] += pData->weight_ptr()[iObs]*( 1 - (2*pData->y_ptr()[iObs]-1)*adF[iObs])*( 1 - (2*pData->y_ptr()[iObs]-1)*adF[iObs]);
	  }
        } // close if(afInBag[iObs
    }
  
  for(iNode=0; iNode<cTermNodes; iNode++)
    {
      if(pTreeComps->GetTermNodes()[iNode]!=NULL)
        {
	  if(vecdDen[iNode] == 0)
            {
	      pTreeComps->GetTermNodes()[iNode]->dPrediction = 0.0;
            }
	  else
            {
	      pTreeComps->GetTermNodes()[iNode]->dPrediction =
		vecdNum[iNode]/vecdDen[iNode];
            }
        }
    }
}


double CHuberized::BagImprovement
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

            if( (2*data.y_ptr()[i]-1)*dF < -1 ){
               dReturnValue += data.weight_ptr()[i]*
                   (-4*(2*data.y_ptr()[i]-1)*dF -
                    -4*(2*data.y_ptr()[i]-1)*(dF+shrinkage*adFadj[i]));
               dW += data.weight_ptr()[i];
            }
            else if ( 1 - (2*data.y_ptr()[i]-1)*dF < 0 ){
               dReturnValue += 0;
               dW += data.weight_ptr()[i];
            }
            else {
               dReturnValue += data.weight_ptr()[i] *
                  ( ( 1 - (2*data.y_ptr()[i]-1)*dF )*( 1 - (2*data.y_ptr()[i]-1)*dF ) -
                    ( 1 - (2*data.y_ptr()[i]-1)*(dF+shrinkage*adFadj[i]) )*( 1 - (2*data.y_ptr()[i]-1)*(dF+shrinkage*adFadj[i]) )
                  );
            }
        }
    }

    return dReturnValue/dW;
}
