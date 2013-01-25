// GBM by Greg Ridgeway  Copyright (C) 2003
// huberized.ccp & huberized.h and associated R code added
// by Harry Southworth, April 2009.

#include "huberized.h"

CHuberized::CHuberized()
{
}

CHuberized::~CHuberized()
{
}


GBMRESULT CHuberized::ComputeWorkingResponse
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adF, 
    double *adZ, 
    double *adWeight,
    bool *afInBag,
    unsigned long nTrain,
    int cIdxOff
)
{
   unsigned long i = 0;
   double dF = 0.0;
   
   for(i=0; i<nTrain; i++)
   {
      dF = adF[i] + ((adOffset==NULL) ? 0.0 : adOffset[i]);
      if( (2*adY[i]-1)*dF < -1)
      {
         adZ[i] = -4 * (2*adY[i]-1);
      }
      else if ( 1 - (2*adY[i]-1)*dF < 0 )
      {
         adZ[i] = 0;
      }
      else
      {
         adZ[i] = -2 * (2*adY[i]-1) * ( 1 - (2*adY[i]-1)*dF );
      }
   }
   return GBM_OK;
}

GBMRESULT CHuberized::InitF
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adWeight,
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
        if(adY[i]==1.0)
        {
            dNum += adWeight[i];
        }
        else
        {
            dDen += adWeight[i];
        }
    }

    dInitF = dNum/dDen;

    return GBM_OK;
}


double CHuberized::Deviance
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adWeight,
    double *adF,
    unsigned long cLength,
    int cIdxOff
)
{
   unsigned long i=0;
   double dL = 0.0;
   double dF = 0.0;
   double dW = 0.0;

   if(adOffset==NULL)
   {
      for(i=cIdxOff; i<cLength+cIdxOff; i++)
      {
        if ( (2*adY[i]-1)*adF[i] < -1 )
         {
            dL += -adWeight[i]*4*(2*adY[i]-1)*adF[i];
            dW += adWeight[i];
         }
         else if ( 1 - (2*adY[i]-1)*adF[i] < 0 ){
            dL += 0;
            dW += adWeight[i];
         }
         else {
            dL += adWeight[i]*( 1 - (2*adY[i]-1)*adF[i] )*( 1 - (2*adY[i]-1)*adF[i] );
            dW += adWeight[i];
         }
      }
   } // close if (adOffset==NULL)
   else
   {
      for(i=cIdxOff; i<cLength+cIdxOff; i++)
      {
         dF = adOffset[i]+adF[i];
         if ( (2*adY[i]-1)*adF[i] < -1 )
         {
            dL += -adWeight[i]*4*(2*adY[i]-1)*dF;
            dW += adWeight[i];
         }
         else if ( 1 - (2*adY[i]-1)*dF < 0 )
         {
            dL += 0;
            dW += adWeight[i];
         }
         else
         {
            dL += adWeight[i] * ( 1 - (2*adY[i]-1)*dF ) * 
                                ( 1 - (2*adY[i]-1)*dF );
            dW += adWeight[i];
         }
      } // close for(
   } // close else

   return dL/dW;
}


GBMRESULT CHuberized::FitBestConstant
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adW,
    double *adF,
    double *adZ,
    unsigned long *aiNodeAssign,
    unsigned long nTrain,
    VEC_P_NODETERMINAL vecpTermNodes,
    unsigned long cTermNodes,
    unsigned long cMinObsInNode,
    bool *afInBag,
    double *adFadj,
   int cIdxOff
)
{
    GBMRESULT hr = GBM_OK;

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
           dF = adF[iObs] + ((adOffset==NULL) ? 0.0 : adOffset[iObs]);
           if( (2*adY[iObs]-1)*adF[iObs] < -1 ){ 
              vecdNum[aiNodeAssign[iObs]] +=
                adW[iObs]*4*(2*adY[iObs]-1);
              vecdDen[aiNodeAssign[iObs]] += 
                -adW[iObs]*4*(2*adY[iObs]-1)*dF;
           }
           else if ( 1 - (2*adY[iObs]-1)*adF[iObs] < 0 ){
              vecdNum[aiNodeAssign[iObs]] += 0;
              vecdDen[aiNodeAssign[iObs]] += 0;
           }
           else{
              vecdNum[aiNodeAssign[iObs]] += adW[iObs]*2*(2*adY[iObs]-1)*( 1 - (2*adY[iObs]-1)*adF[iObs] );
              vecdDen[aiNodeAssign[iObs]] += adW[iObs]*( 1 - (2*adY[iObs]-1)*adF[iObs])*( 1 - (2*adY[iObs]-1)*adF[iObs]);
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

    return hr;
}


double CHuberized::BagImprovement
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adWeight,
    double *adF,
    double *adFadj,
    bool *afInBag,
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
            dF = adF[i] + ((adOffset==NULL) ? 0.0 : adOffset[i]);

            if( (2*adY[i]-1)*dF < -1 ){
               dReturnValue += adWeight[i]*
                   (-4*(2*adY[i]-1)*dF -
                    -4*(2*adY[i]-1)*(dF+dStepSize*adFadj[i]));
               dW += adWeight[i];
            }
            else if ( 1 - (2*adY[i]-1)*dF < 0 ){
               dReturnValue += 0;
               dW += adWeight[i];
            }
            else {
               dReturnValue += adWeight[i] *
                  ( ( 1 - (2*adY[i]-1)*dF )*( 1 - (2*adY[i]-1)*dF ) -
                    ( 1 - (2*adY[i]-1)*(dF+dStepSize*adFadj[i]) )*( 1 - (2*adY[i]-1)*(dF+dStepSize*adFadj[i]) )
                  );
            }
        }
    }

    return dReturnValue/dW;
}
