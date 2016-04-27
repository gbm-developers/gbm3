//------------------------------------------------------------------------------
//
//  File:       censoredCoxState.h
//
//  Description: censored CoxPH methods
//
//	Author: 	James Hickey
//------------------------------------------------------------------------------

#ifndef __censoredCoxState_h__
#define __censoredCoxState_h__

//------------------------------
// Includes
//------------------------------
#include "dataset.h"
#include "distribution.h"
#include "genericCoxState.h"
#include <Rcpp.h>

//------------------------------
// Class Definition
//------------------------------
class CensoredCoxState: public GenericCoxState
{
public:
	//----------------------
	// Public Constructors
	//----------------------
	CensoredCoxState(CCoxPH* coxPhPtr): coxPh(coxPhPtr){};

	//---------------------
	// Public destructor
	//---------------------
	~CensoredCoxState(){coxPh = NULL;};

	//---------------------
	// Public Functions
	//---------------------
	void ComputeWorkingResponse
	(
		const CDataset* pData,
	    const double *adF,
	    double *adZ
	)
	{
		// Initialize parameters
		std::vector<double> martingaleResid(pData->get_trainSize(), 0.0);
		double loglik = LogLikelihood(pData->get_trainSize(), pData, adF, &martingaleResid[0], false);

		// Fill up response
		for(long i = 0; i < pData->get_trainSize(); i++)
		{
			if(pData->GetBagElem(i))
			{
				adZ[i] = martingaleResid[i]; // From chain rule
			}
		}
	}

	void FitBestConstant
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

		vector<double> vecdP;
		vector<double> vecdG;
		vector<unsigned long> veciK2Node(cTermNodes, 0);
		vector<unsigned long> veciNode2K(cTermNodes, 0);

		matrix<double> matH;
		matrix<double> matHinv;

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

				if(coxPh->StatusVec()[i]==1.0)
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

	double Deviance
	(
		const long cLength,
		const CDataset* pData,
	    const double *adF
	)
	{
		// Initialize Parameters
	    double loglik = 0.0;
	    std::vector<double> martingaleResid(cLength, 0.0);

	    // Calculate Deviance
	    loglik = LogLikelihood(cLength, pData, adF, &martingaleResid[0]);

	    return -loglik;
	}

	double BagImprovement
	(
		const CDataset& data,
		const double *adF,
		const bag& afInBag,
	  const double shrinkage,
	  const double* adFadj
	)
	{
		// Initialize Parameters
	    double loglikeNoAdj = 0.0;
	    double loglikeWithAdj = 0.0;
	    std::vector<double> martingaleResidNoAdj(data.get_trainSize(), 0.0);
	    std::vector<double> martingaleResidWithAdj(data.get_trainSize(), 0.0);
	    std::vector<double> etaAdj(data.get_trainSize(), 0.0);

	    // Fill up the adjusted and shrunk eta
	    for(long i = 0; i < data.get_trainSize(); i++)
	    {
	    	if(!data.GetBagElem(i))
	    	{
		    	etaAdj[i] = adF[i] + shrinkage * adFadj[i];

	    	}
	    	else
	    	{
	    		etaAdj[i]  = adF[i];
	    	}
	    }

	    // Calculate likelihoods - data not in bags
	    loglikeNoAdj = LogLikelihood(data.get_trainSize(), &data, adF, &martingaleResidNoAdj[0], false, false);
	    loglikeWithAdj = LogLikelihood(data.get_trainSize(), &data, &etaAdj[0], &martingaleResidWithAdj[0], false, false);

	    return (loglikeWithAdj - loglikeNoAdj);
	}

private:
	CCoxPH* coxPh;

	double LogLikelihood(const int n, const CDataset* pData, const double* eta,
								double* resid, bool skipBag=true, bool checkInBag=true)
	{
	    int i, j, k, ksave;
	    int person, p2;
	    int istrat, indx1;   /* this counts up over the strata */
	    double cumhaz, hazard;
	    int nrisk, ndeath;
	    double deathwt, denom, temp, center;
	    double esum, dtime, e_hazard;
	    double loglik, d_denom;

	    /*
	    **  'person' walks through the the data from 1 to n, p2= sort2[person].
	    **     sort2[0] points to the largest stop time, sort2[1] the next, ...
	    **  'dtime' is a scratch variable holding the time of current interest
	    */
	    istrat=0; indx1=0;
	    denom =0;  /* S in the math explanation */
	    cumhaz =0;
	    nrisk =0;   /* number at risk */
	    esum =0;  /*cumulative eta, used for rescaling */
	    loglik =0;
	    center = eta[coxPh->EndTimeIndices()[0]] + pData->offset_ptr(false)[coxPh->EndTimeIndices()[0]];
		//throw GBM::failure("Test CoxPH");

	    for (person=0; person<n; )
	    {

	    	// Check if bagging is required
	    	if(skipBag || (pData->GetBagElem(person)==checkInBag))
	    	{
	    		p2 = coxPh->EndTimeIndices()[person];
	    		if (coxPh->StatusVec()[p2] ==0)
				{
					/* add the subject to the risk set */
					resid[p2] = exp(eta[p2] + pData->offset_ptr(false)[p2] - center) * cumhaz;
					nrisk++;
					denom  += pData->weight_ptr()[p2]* exp(eta[p2] + pData->offset_ptr(false)[p2] - center);
					esum += eta[p2] + pData->offset_ptr(false)[p2];
					person++;
				}
				else
				{
					dtime = pData->y_ptr()[p2];  /* found a new, unique death time */
					/*
					**        Add up over this death time, for all subjects
					*/
					ndeath =0;
					deathwt =0;
					d_denom =0;  /*contribution to denominator by death at dtime */
					for (k=person; k<coxPh->StrataVec()[istrat]; k++)
					{
						// Check in loop over stratum that person has correct bag
						// properties
						if(skipBag || (pData->GetBagElem(k)==checkInBag))
						{
							p2 = coxPh->EndTimeIndices()[k];
							if (pData->y_ptr()[p2]  < dtime) break;  /* only tied times */

							nrisk++;
							denom += pData->weight_ptr()[p2] * exp(eta[p2] + pData->offset_ptr(false)[p2] - center);
							esum += eta[p2] + pData->offset_ptr(false)[p2];
							if (coxPh->StatusVec()[p2] ==1)
							{
								ndeath ++;
								deathwt += pData->weight_ptr()[p2];
								d_denom += pData->weight_ptr()[p2] * exp(eta[p2] + pData->offset_ptr(false)[p2] - center);
								loglik  += pData->weight_ptr()[p2]*(eta[p2] + pData->offset_ptr(false)[p2] - center);
							 }
						}
					}
					ksave = k;
					/* compute the increment in hazard
					** hazard = usual increment
					** e_hazard = efron increment, for tied deaths only
					*/
					if (coxPh->TieApproxMethod()==0 || ndeath==1)
					{ /* Breslow */
						loglik -= deathwt* log(denom);
						hazard = deathwt /denom;
						e_hazard = hazard;
					}
					else
					{ /* Efron */
						hazard =0;
						e_hazard =0;  /* hazard experienced by a tied death */
						deathwt /= ndeath;   /* average weight of each death */
						for (k=0; k <ndeath; k++)
						{
							if(skipBag || (pData->GetBagElem(k) == checkInBag))
							{
								temp = (double)k /ndeath;    /* don't do integer division*/
								loglik -= deathwt * log(denom - temp*d_denom);
								hazard += deathwt/(denom - temp*d_denom);
								e_hazard += (1-temp) *deathwt/(denom - temp*d_denom);
							}

						}
					}

					/* Give initial value to all intervals ending at this time
					** If tied censors are sorted before deaths (which at least some
					**  callers of this routine do), then the else below will never
					**  occur.
					*/
					temp = cumhaz + (hazard -e_hazard);
					for (; person < ksave; person++)
					{
						// Check if person in/out of bag
						if(skipBag || (pData->GetBagElem(person)==checkInBag))
						{
							p2 = coxPh->EndTimeIndices()[person];
							if (coxPh->StatusVec()[p2] ==1) resid[p2] = 1 + temp*exp(eta[p2] + pData->offset_ptr(false)[p2] - center);
							else resid[p2] = cumhaz * exp(eta[p2] + pData->offset_ptr(false)[p2] - center);
						}

					}
					cumhaz += hazard;

					/* see if we need to shift the centering (a rare case) */
					if (fabs(esum/nrisk - center) > recenter)
					{
						temp = esum/nrisk - center;
						center += temp;
						denom /=  exp(temp);
					}
				}

				/* clean up at the end of a strata */
				if (person == coxPh->StrataVec()[istrat])
				{
					for (; indx1< coxPh->StrataVec()[istrat]; indx1++)
					{
						// Check bagging status
						if(skipBag || (pData->GetBagElem(indx1)==checkInBag))
						{
							p2 = coxPh->EndTimeIndices()[indx1];
							resid[p2] -= cumhaz * exp(eta[p2] + pData->offset_ptr(false)[p2] - center);
						}
					}
					cumhaz =0;
					denom = 0;
					istrat++;
				}
	    	}
	    	else
			{
				// Increment person if not in bag
				person++;
			}

	    }
	    return(loglik);
	}

};
#endif //__censoredCoxState_h__
