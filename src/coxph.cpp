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
									int& cGroups, int& cTrain)
{
	return new CCoxPH(radMisc);
}

CCoxPH::~CCoxPH()
{
}


void CCoxPH::ComputeWorkingResponse
////////////////////////////////////////////
// weightedQuantile
//
// Function to return the weighted quantile of
// a vector of a given length
//
// Parameters: iN     - Length of vector
//             adV    - Vector of doubles
//             pData->weight_ptr()    - Array of weights
//             dAlpha - Quantile to calculate (0.5 for median)
//
// Returns :   Weighted quantile
/////////////////////////////////////////////////ponse
(
	const CDataset* pData,
    const double *adF,
    double *adZ,
    const bag& afInBag,
    unsigned long nTrain
)
{
    unsigned long i = 0;
    double dF = 0.0;
    double dTot = 0.0;
    double dRiskTot = 0.0;

    vecdRiskTot.resize(nTrain);
    dRiskTot = 0.0;
    for(i=0; i<nTrain; i++)
    {
        if(afInBag[i])
        {
            dF = adF[i] + ((pData->offset_ptr(false)==NULL) ? 0.0 : pData->offset_ptr(false)[i]);
            dRiskTot += pData->weight_ptr()[i]*std::exp(dF);
            vecdRiskTot[i] = dRiskTot;
        }
    }

    dTot = 0.0;
    for(i=nTrain-1; i!=ULONG_MAX; i--) // i is unsigned so wraps to ULONG_MAX
    {
        if(afInBag[i])
        {
            if(adDelta[i]==1.0)
            {
                dTot += pData->weight_ptr()[i]/vecdRiskTot[i];
            }
            dF = adF[i] + ((pData->offset_ptr(false)==NULL) ? 0.0 : pData->offset_ptr(false)[i]);
            adZ[i] = adDelta[i] - std::exp(dF)*dTot;
        }
    }

}



void CCoxPH::InitF
(
	const CDataset* pData,
    double &dInitF,
    unsigned long cLength
)
{
    dInitF = 0.0;
}


double CCoxPH::Deviance
(
	const CDataset* pData,
    const double *adF,
    unsigned long cLength,
    bool isValidationSet
)
{
    unsigned long i=0;
    double dL = 0.0;
    double dF = 0.0;
    double dW = 0.0;
    double dTotalAtRisk = 0.0;

    if(isValidationSet)
    {
    	pData->shift_to_validation();
    	adDelta = shift_ptr(CDistribution::misc_ptr(false), pData->get_trainSize());
    }

    dTotalAtRisk = 0.0; 
    for(i=0; i!=cLength; i++)
    {
        dF = adF[i] + ((pData->offset_ptr(false)==NULL) ? 0.0 : pData->offset_ptr(false)[i]);
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
    return -2*dL/dW;
}


void CCoxPH::FitBestConstant
(
	const CDataset* pData,
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
        if(vecpTermNodes[i]->cN >= cMinObsInNode)
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
    for(i=0; i<nTrain; i++)
    {
        if(afInBag[i] && (vecpTermNodes[aiNodeAssign[i]]->cN >= cMinObsInNode))
        {
            dF = adF[i] + ((pData->offset_ptr(false)==NULL) ? 0.0 : pData->offset_ptr(false)[i]);
            vecdP[veciNode2K[aiNodeAssign[i]]] += pData->weight_ptr()[i]*std::exp(dF);
            dRiskTot += pData->weight_ptr()[i]*std::exp(dF);

            if(adDelta[i]==1.0)
            {
                // compute g and H
                for(k=0; k<K-1; k++)
                {
                    vecdG[k] +=
                        pData->weight_ptr()[i]*((aiNodeAssign[i]==veciK2Node[k]) - vecdP[k]/dRiskTot);

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
        vecpTermNodes[k]->dPrediction = 0.0;
    }
    for(m=0; m<K-1; m++)
    {
        for(k=0; k<K-1; k++)
        {
            matH.getvalue(k,m,dTemp,fTemp);
            if(!R_FINITE(dTemp)) // occurs if matH was not invertible
            {
                vecpTermNodes[veciK2Node[k]]->dPrediction = 0.0;
                break;
            }
            else
            {
                vecpTermNodes[veciK2Node[k]]->dPrediction -= dTemp*vecdG[m];
            }
          }
    }
    // vecpTermNodes[veciK2Node[K-1]]->dPrediction = 0.0; // already set to 0.0
}


double CCoxPH::BagImprovement
(
	const CDataset* pData,
    const double *adF,
    const double *adFadj,
    const bag& afInBag,
    double dStepSize,
    unsigned long nTrain
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
    for(i=0; i<nTrain; i++)
    {
        if(!afInBag[i])
        {
            dNum += pData->weight_ptr()[i]*std::exp(dF + dStepSize*adFadj[i]);
            dDen += pData->weight_ptr()[i]*std::exp(dF);
            if(adDelta[i]==1.0)
            {
                dReturnValue +=
                    pData->weight_ptr()[i]*(dStepSize*adFadj[i] - std::log(dNum) + log(dDen));
                dW += pData->weight_ptr()[i];
            }
        }
    }

    return dReturnValue/dW;
}



