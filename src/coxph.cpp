//-----------------------------------
//
// File: coxph.cpp
//
// Description: Cox partial hazards model.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "coxph.h"

//----------------------------------------
// Function Members - Private
//----------------------------------------
CCoxPH::CCoxPH(SEXP radMisc): CDistribution(radMisc)
{
	adDelta = CDistribution::misc_ptr(false);
}


//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CCoxPH::Create(SEXP radMisc,
							  const char* szIRMeasure,
							  int& cTrain)
{
	return new CCoxPH(radMisc);
}

CCoxPH::~CCoxPH()
{
}


void CCoxPH::ComputeWorkingResponse
(
	const CDataset* pData,
    const double *adF,
    double *adZ
)
{
    double dF = 0.0;
    double dTot = 0.0;
    double dRiskTot = 0.0;

    vecdRiskTot.resize(pData->get_trainSize());
    dRiskTot = 0.0;
    for(unsigned long i=0; i<pData->get_trainSize(); i++)
    {
        if(pData->GetBagElem(i))

        {
            dF = adF[i] +  pData->offset_ptr(false)[i];
            dRiskTot += pData->weight_ptr()[i]*std::exp(dF);
            vecdRiskTot[i] = dRiskTot;
        }
    }
    dTot = 0.0;
    for(long i= pData->get_trainSize()-1; i != -1; i--)
    {
        if(pData->GetBagElem(i))
        {
            if(adDelta[i]==1.0)
            {
                dTot += pData->weight_ptr()[i]/vecdRiskTot[i];
            }
            dF = adF[i] +  pData->offset_ptr(false)[i];
            adZ[i] = adDelta[i] - std::exp(dF)*dTot;
        }
    }
}



double CCoxPH::InitF
(
	const CDataset* pData
)
{
    return 0.0;
}


double CCoxPH::Deviance
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
    double dTotalAtRisk = 0.0;

    long cLength = pData->get_trainSize();
    if(isValidationSet)
    {
    	pData->shift_to_validation();
    	adDelta = shift_ptr(CDistribution::misc_ptr(false), pData->get_trainSize());
    	cLength = pData->GetValidSize();
    }

    dTotalAtRisk = 0.0; 
    for(i=0; i!=cLength; i++)
    {
        dF = adF[i] +  pData->offset_ptr(false)[i];
        dTotalAtRisk += pData->weight_ptr()[i]*std::exp(dF);
        if(adDelta[i]==1.0)
        {
            dL += pData->weight_ptr()[i]*(dF - std::log(dTotalAtRisk));
            dW += pData->weight_ptr()[i];
        }
    }

    // Shift back for future calculations if required
    if(isValidationSet)
    {
    	pData->shift_to_train();
    	adDelta = CDistribution::misc_ptr(false);
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


void CCoxPH::FitBestConstant
(
	const CDataset* pData,
    const double *adF,
    unsigned long cTermNodes,
    double* adZ,
    CTreeComps* pTreeComps
)
{
    double dF = 0.0;
    double dRiskTot = 0.0;
    unsigned long i = 0;
    unsigned long k = 0;
    unsigned long m = 0;

    double dTemp = 0.0;
    bool fTemp = false;
    unsigned long K = 0;
    veciK2Node.resize(cTermNodes);
    veciNode2K.resize(cTermNodes);

    for(i=0; i<cTermNodes; i++)
    {
        veciNode2K[i] = 0;
        if(pTreeComps->GetTermNodes()[i]->cN >= pTreeComps->GetMinNodeObs())
        {
            veciK2Node[K] = i;
            veciNode2K[i] = K;
            K++;
        }
    }

    vecdP.resize(K);

    matH.setactualsize(K-1);
    vecdG.resize(K-1);
    vecdG.assign(K-1,0.0);

    // zero the Hessian
    for(k=0; k<K-1; k++)
    {
        for(m=0; m<K-1; m++)
        {
            matH.setvalue(k,m,0.0);
        }
    }

    // get the gradient & Hessian, Ridgeway (1999) pp. 100-101
    // correction from Ridgeway (1999): fix terminal node K-1 prediction to 0.0
    //      for identifiability
    dRiskTot = 0.0;
    vecdP.assign(K,0.0);
    for(i=0; i<pData->get_trainSize(); i++)
    {
        if(pData->GetBagElem(i) && (pTreeComps->GetTermNodes()[pTreeComps->GetNodeAssign()[i]]->cN >= pTreeComps->GetMinNodeObs()))
        {
            dF = adF[i] + ((pData->offset_ptr(false)==NULL) ? 0.0 : pData->offset_ptr(false)[i]);
            vecdP[veciNode2K[pTreeComps->GetNodeAssign()[i]]] += pData->weight_ptr()[i]*std::exp(dF);
            dRiskTot += pData->weight_ptr()[i]*std::exp(dF);

            if(adDelta[i]==1.0)
            {
                // compute g and H
                for(k=0; k<K-1; k++)
                {
                    vecdG[k] +=
                        pData->weight_ptr()[i]*((pTreeComps->GetNodeAssign()[i]==veciK2Node[k]) - vecdP[k]/dRiskTot);

                    matH.getvalue(k,k,dTemp,fTemp);
                    matH.setvalue(k,k,dTemp -
                        pData->weight_ptr()[i]*vecdP[k]/dRiskTot*(1-vecdP[k]/dRiskTot));
                    for(m=0; m<k; m++)
                    {
                        matH.getvalue(k,m,dTemp,fTemp);
                        dTemp += pData->weight_ptr()[i]*vecdP[k]/dRiskTot*vecdP[m]/dRiskTot;
                        matH.setvalue(k,m,dTemp);
                        matH.setvalue(m,k,dTemp);
                    }
                }
            }
        }
    }

    /*
    for(k=0; k<K-1; k++)
    {
        for(m=0; m<K-1; m++)
        {
            matH.getvalue(k,m,dTemp,fTemp);
            Rprintf("%f ",dTemp);
        }
        Rprintf("\n");
    }
    */

    // one step to get leaf predictions
    matH.invert();

    for(k=0; k<cTermNodes; k++)
    {
        pTreeComps->GetTermNodes()[k]->dPrediction = 0.0;
    }
    for(m=0; m<K-1; m++)
    {
        for(k=0; k<K-1; k++)
        {
            matH.getvalue(k,m,dTemp,fTemp);
            if(!R_FINITE(dTemp)) // occurs if matH was not invertible
            {
                pTreeComps->GetTermNodes()[veciK2Node[k]]->dPrediction = 0.0;
                break;
            }
            else
            {
                pTreeComps->GetTermNodes()[veciK2Node[k]]->dPrediction -= dTemp*vecdG[m];
            }
          }
    }
    // vecpTermNodes[veciK2Node[K-1]]->dPrediction = 0.0; // already set to 0.0
}


double CCoxPH::BagImprovement
(
	const CDataset& data,
    const double *adF,
    const bag& afInBag,
  const double shrinkage,
  const double* adFadj
)
{
    double dReturnValue = 0.0;
    double dNum = 0.0;
    double dDen = 0.0;
    double dF = 0.0;
    double dW = 0.0;
    unsigned long i = 0;

    dNum = 0.0;
    dDen = 0.0;
    for(i=0; i< data.get_trainSize(); i++)
    {
        if(!data.GetBagElem(i))
        {
            dNum += data.weight_ptr()[i]*std::exp(dF + shrinkage*adFadj[i]);
            dDen += data.weight_ptr()[i]*std::exp(dF);
            if(adDelta[i]==1.0)
            {
                dReturnValue +=
                    data.weight_ptr()[i]*(shrinkage*adFadj[i] - std::log(dNum) + log(dDen));
                dW += data.weight_ptr()[i];
            }
        }
    }

    //std::cout << dReturnValue/dW << endl;
    return dReturnValue/dW;
}



