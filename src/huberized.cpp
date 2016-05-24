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
CHuberized::CHuberized()
{
}

//----------------------------------------
// Function Members - Public
//----------------------------------------

CDistribution* CHuberized::Create(DataDistParams& distParams)
{
	return new CHuberized();
}

CHuberized::~CHuberized()
{
}


void CHuberized::ComputeWorkingResponse
(
	const CDataset& data,
    const double *adF,
    double *adZ
)
{
   unsigned long i = 0;
   double dF = 0.0;

   for(i=0; i<data.get_trainsize(); i++)
   {
      dF = adF[i] + data.offset_ptr()[i];
      if( (2*data.y_ptr()[i]-1)*dF < -1)
      {
         adZ[i] = -4 * (2*data.y_ptr()[i]-1);
      }
      else if ( 1 - (2*data.y_ptr()[i]-1)*dF < 0 )
      {
         adZ[i] = 0;
      }
      else
      {
         adZ[i] = -2 * (2*data.y_ptr()[i]-1) * ( 1 - (2*data.y_ptr()[i]-1)*dF );
      }
   }
}

double CHuberized::InitF
(
	const CDataset& data
)
{
    unsigned long i=0;
    double dNum = 0.0;
    double dDen = 0.0;

    for(i=0; i<data.get_trainsize(); i++)
    {
        if(data.y_ptr()[i]==1.0)
        {
            dNum += data.weight_ptr()[i];
        }
        else
        {
            dDen += data.weight_ptr()[i];
        }
    }

    return dNum/dDen;
}


double CHuberized::Deviance
(
	const CDataset& data,
    const double *adF,
    bool isValidationSet
)
{
   unsigned long i=0;
   double dL = 0.0;
   double dF = 0.0;
   double dW = 0.0;

   unsigned long cLength = data.get_trainsize();
   if(isValidationSet)
   {
	   data.shift_to_validation();
	   cLength = data.get_validsize();
   }


  for(i=0; i<cLength; i++)
  {
	 dF = data.offset_ptr()[i]+adF[i];
	 if ( (2*data.y_ptr()[i]-1)*adF[i] < -1 )
	 {
		dL += -data.weight_ptr()[i]*4*(2*data.y_ptr()[i]-1)*dF;
		dW += data.weight_ptr()[i];
	 }
	 else if ( 1 - (2*data.y_ptr()[i]-1)*dF < 0 )
	 {
		dL += 0;
		dW += data.weight_ptr()[i];
	 }
	 else
	 {
		dL += data.weight_ptr()[i] * ( 1 - (2*data.y_ptr()[i]-1)*dF ) *
							( 1 - (2*data.y_ptr()[i]-1)*dF );
		dW += data.weight_ptr()[i];
	 }
  } // close for(


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


void CHuberized::FitBestConstant
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
  
  vector<double> vecdNum(cTermNodes, 0.0);
  vector<double> vecdDen(cTermNodes, 0.0);

  for(iObs=0; iObs<data.get_trainsize(); iObs++)
    {
      if(data.get_bag_element(iObs))
        {
	  dF = adF[iObs] +  data.offset_ptr()[iObs];
	  if( (2*data.y_ptr()[iObs]-1)*adF[iObs] < -1 )
	  {
	    vecdNum[treeComps.get_node_assignments()[iObs]] +=
	      data.weight_ptr()[iObs]*4*(2*data.y_ptr()[iObs]-1);
	    vecdDen[treeComps.get_node_assignments()[iObs]] +=
	      -data.weight_ptr()[iObs]*4*(2*data.y_ptr()[iObs]-1)*dF;
	  }
	  else if ( 1 - (2*data.y_ptr()[iObs]-1)*adF[iObs] < 0 ){
	    vecdNum[treeComps.get_node_assignments()[iObs]] += 0;
	    vecdDen[treeComps.get_node_assignments()[iObs]] += 0;
	  }
	  else{
	    vecdNum[treeComps.get_node_assignments()[iObs]] += data.weight_ptr()[iObs]*2*(2*data.y_ptr()[iObs]-1)*( 1 - (2*data.y_ptr()[iObs]-1)*adF[iObs] );
	    vecdDen[treeComps.get_node_assignments()[iObs]] += data.weight_ptr()[iObs]*( 1 - (2*data.y_ptr()[iObs]-1)*adF[iObs])*( 1 - (2*data.y_ptr()[iObs]-1)*adF[iObs]);
	  }
        } // close if(afInBag[iObs
    }
  
  for(iNode=0; iNode<cTermNodes; iNode++)
    {
      if(treeComps.get_terminal_nodes()[iNode]!=NULL)
        {
	  if(vecdDen[iNode] == 0)
            {
	      treeComps.get_terminal_nodes()[iNode]->prediction = 0.0;
            }
	  else
            {
	      treeComps.get_terminal_nodes()[iNode]->prediction =
		vecdNum[iNode]/vecdDen[iNode];
            }
        }
    }
}


double CHuberized::BagImprovement
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
            dF = adF[i] +  data.offset_ptr()[i];

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
               //TODO: Does this require an dW+= ?
            }
        }
    }

    //TODO: Check if weights are all zero for validation set
   if((dW == 0.0) && (dReturnValue == 0.0))
   {
	   return nan("");
   }
   else if(dW == 0.0)
   {

	   return copysign(HUGE_VAL, dReturnValue);
   }

    return dReturnValue/dW;
}
