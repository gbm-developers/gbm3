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
CHuberized::CHuberized(SEXP radMisc, const CDataset& data): CDistribution(radMisc, data)
{
}

//----------------------------------------
// Function Members - Public
//----------------------------------------

std::auto_ptr<CDistribution> CHuberized::Create(SEXP radMisc, const CDataset& data,
										const char* szIRMeasure,
										int& cGroups, int& cTrain)
{
	return std::auto_ptr<CDistribution>(new CHuberized(radMisc, data));
}

CHuberized::~CHuberized()
{
}


void CHuberized::ComputeWorkingResponse
(
    const double *adF,
    double *adZ,
    const bag& afInBag,
    unsigned long nTrain
)
{
   unsigned long i = 0;
   double dF = 0.0;

   for(i=0; i<nTrain; i++)
   {
      dF = adF[i] + ((pData->offset_ptr(false)==NULL) ? 0.0 : pData->offset_ptr(false)[i]);
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

void CHuberized::InitF
(
    double &dInitF,
    unsigned long cLength
)
{
    unsigned long i=0;
    double dNum = 0.0;
    double dDen = 0.0;

    dInitF = 0.0;

    for(i=0; i<cLength; i++)
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

    dInitF = dNum/dDen;
}


double CHuberized::Deviance
(
    const double *adF,
    unsigned long cLength,
    bool isValidationSet
)
{
   unsigned long i=0;
   double dL = 0.0;
   double dF = 0.0;
   double dW = 0.0;

   if(isValidationSet)
   {
	   pData->shift_to_validation();
   }

   if(pData->offset_ptr(false)==NULL)
   {
      for(i=0; i<cLength; i++)
      {
        if ( (2*pData->y_ptr()[i]-1)*adF[i] < -1 )
         {
            dL += -pData->weight_ptr()[i]*4*(2*pData->y_ptr()[i]-1)*adF[i];
            dW += pData->weight_ptr()[i];
         }
         else if ( 1 - (2*pData->y_ptr()[i]-1)*adF[i] < 0 ){
            dL += 0;
            dW += pData->weight_ptr()[i];
         }
         else {
            dL += pData->weight_ptr()[i]*( 1 - (2*pData->y_ptr()[i]-1)*adF[i] )*( 1 - (2*pData->y_ptr()[i]-1)*adF[i] );
            dW += pData->weight_ptr()[i];
         }
      }
   } // close if (pData->offset_ptr(false)==NULL)
   else
   {
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
   } // close else

   if(isValidationSet)
   {
	   pData->shift_to_train();
   }

   return dL/dW;
}


void CHuberized::FitBestConstant
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
  vecdNum.resize(cTermNodes);
  vecdNum.assign(vecdNum.size(),0.0);
  vecdDen.resize(cTermNodes);
  vecdDen.assign(vecdDen.size(),0.0);
  
  for(iObs=0; iObs<nTrain; iObs++)
    {
      if(afInBag[iObs])
        {
	  dF = adF[iObs] + ((pData->offset_ptr(false)==NULL) ? 0.0 : pData->offset_ptr(false)[iObs]);
	  if( (2*pData->y_ptr()[iObs]-1)*adF[iObs] < -1 ){
	    vecdNum[aiNodeAssign[iObs]] +=
	      pData->weight_ptr()[iObs]*4*(2*pData->y_ptr()[iObs]-1);
	    vecdDen[aiNodeAssign[iObs]] +=
	      -pData->weight_ptr()[iObs]*4*(2*pData->y_ptr()[iObs]-1)*dF;
	  }
	  else if ( 1 - (2*pData->y_ptr()[iObs]-1)*adF[iObs] < 0 ){
	    vecdNum[aiNodeAssign[iObs]] += 0;
	    vecdDen[aiNodeAssign[iObs]] += 0;
	  }
	  else{
	    vecdNum[aiNodeAssign[iObs]] += pData->weight_ptr()[iObs]*2*(2*pData->y_ptr()[iObs]-1)*( 1 - (2*pData->y_ptr()[iObs]-1)*adF[iObs] );
	    vecdDen[aiNodeAssign[iObs]] += pData->weight_ptr()[iObs]*( 1 - (2*pData->y_ptr()[iObs]-1)*adF[iObs])*( 1 - (2*pData->y_ptr()[iObs]-1)*adF[iObs]);
	  }
        } // close if(afInBag[iObs
    }
  
  for(iNode=0; iNode<cTermNodes; iNode++)
    {
      if(vecpTermNodes[iNode]!=NULL)
        {
	  if(vecdDen[iNode] == 0)
            {
	      vecpTermNodes[iNode]->dPrediction = 0.0;
            }
	  else
            {
	      vecpTermNodes[iNode]->dPrediction =
		vecdNum[iNode]/vecdDen[iNode];
            }
        }
    }
}


double CHuberized::BagImprovement
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

            if( (2*pData->y_ptr()[i]-1)*dF < -1 ){
               dReturnValue += pData->weight_ptr()[i]*
                   (-4*(2*pData->y_ptr()[i]-1)*dF -
                    -4*(2*pData->y_ptr()[i]-1)*(dF+dStepSize*adFadj[i]));
               dW += pData->weight_ptr()[i];
            }
            else if ( 1 - (2*pData->y_ptr()[i]-1)*dF < 0 ){
               dReturnValue += 0;
               dW += pData->weight_ptr()[i];
            }
            else {
               dReturnValue += pData->weight_ptr()[i] *
                  ( ( 1 - (2*pData->y_ptr()[i]-1)*dF )*( 1 - (2*pData->y_ptr()[i]-1)*dF ) -
                    ( 1 - (2*pData->y_ptr()[i]-1)*(dF+dStepSize*adFadj[i]) )*( 1 - (2*pData->y_ptr()[i]-1)*(dF+dStepSize*adFadj[i]) )
                  );
            }
        }
    }

    return dReturnValue/dW;
}
