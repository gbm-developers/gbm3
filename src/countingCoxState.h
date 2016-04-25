//------------------------------------------------------------------------------
//
//  File:       countingCoxState.h
//
//  Description: counting CoxPH methods
//
//	Author: 	James Hickey
//------------------------------------------------------------------------------

#ifndef __countingCoxState_h__
#define __countingCoxState_h__

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
class CountingCoxState: public GenericCoxState
{
public:
	//----------------------
	// Public Constructors
	//----------------------
	CountingCoxState(CCoxPH* coxPhPtr): coxPh(coxPhPtr){};

	//---------------------
	// Public destructor
	//---------------------
	~CountingCoxState(){coxPh = NULL;};

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
		double loglik = LogLikelihoodTiedTimes(pData->get_trainSize(), pData, adF, &martingaleResid[0], false);

		// Fill up response
		for(long i = 0; i < pData->get_trainSize(); i++)
		{
			if(pData->GetBagElem(i))
			{
				adZ[i] = -martingaleResid[i]; // From chain rule
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
		// Currently does no adjustment
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
		loglik = LogLikelihoodTiedTimes(cLength, pData, adF, &martingaleResid[0]);

		return loglik;
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
			if(data.GetBagElem(i))
			{
				etaAdj[i] = adF[i] + shrinkage * adFadj[i];

			}
		}

		// Calculate likelihoods
		loglikeNoAdj = LogLikelihoodTiedTimes(data.get_trainSize(), &data, adF, &martingaleResidNoAdj[0], false);
		loglikeWithAdj = LogLikelihoodTiedTimes(data.get_trainSize(), &data, &etaAdj[0], &martingaleResidWithAdj[0], false, false);

		return (loglikeWithAdj - loglikeNoAdj);
	}

private:
	CCoxPH* coxPh;
	double LogLikelihoodTiedTimes(const int n, const CDataset* pData, const double* eta,
										  double* resid, bool skipBag=true, bool checkInBag=true)
	{
	    int i, j, k, ksave;
	    int person, p2, indx1, p1;
	    int istrat;
	    double cumhaz, hazard;
	    int nrisk, ndeath;
	    double deathwt, denom, temp, center;
	    double esum, dtime, e_hazard;
	    double loglik, d_denom;
	    int stratastart;    /* the first obs of each stratum */

	    /*
	    **  'person' walks through the the data from 1 to n, p2= sort2[person].
	    **     sort2[0] points to the largest stop time, sort2[1] the next, ...
	    **  'dtime' is a scratch variable holding the time of current interest
	    **  'indx1' walks through the start times.  It will be smaller than
	    **    'person': if person=27 that means that 27 subjects have time2 >=dtime,
	    **    and are thus potential members of the risk set.  If 'indx1' =9,
	    **    that means that 9 subjects have start >=time and thus are NOT part
	    **    of the risk set.  (stop > start for each subject guarrantees that
	    **    the 9 are a subset of the 27). p1 = sort1[indx1]
	    **  Basic algorithm: move 'person' forward, adding the new subject into
	    **    the risk set.  If this is a new, unique death time, take selected
	    **    old obs out of the sums, add in obs tied at this time, then update
	    **    the cumulative hazard. Everything resets at the end of a stratum.
	    **  The sort order is from large time to small, so we encounter a subjects
	    **    ending time first, then their start time.
	    **  The martingale residual for a subject is
	    **     status - (cumhaz at entry - cumhaz at exit)*score
	    */

	    istrat=0;
	    indx1 =0;
	    denom =0;  /* S in the math explanation */
	    cumhaz =0;
	    nrisk =0;   /* number at risk */
	    esum =0;  /*cumulative eta, used for rescaling */
	    center = eta[coxPh->EndTimeIndices()[0]] + pData->offset_ptr(false)[coxPh->EndTimeIndices()[0]];
	    stratastart =0;   /* first strata starts at index 0 */
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
					dtime = pData->y_ptr(1)[p2];  /* found a new, unique death time */

					/*
					** Remove those subjects whose start time is to the right
					**  from the risk set, and finish computation of their residual
					*/
					temp = denom;
					for (;  indx1 <person; indx1++)
					{
						p1 = coxPh->StartTimeIndices()[indx1];
						if (pData->y_ptr(0)[p1] < dtime) break; /* still in the risk set */

						nrisk--;
						resid[p1] -= cumhaz* exp(eta[p1] + pData->offset_ptr(false)[p1] - center);
						denom  -= pData->weight_ptr()[p1] * exp(eta[p1] + pData->offset_ptr(false)[p1] - center);
						esum -= eta[p1] + pData->offset_ptr(false)[p1];
					}
					if (nrisk==0)
					{
						/* everyone was removed!
						** This happens with manufactured start/stop
						**  data sets that have some g(time) as a covariate.
						**  Just like a strata, reset the sums.
						*/
						denom =0;
						esum =0;
					}

					/*
					**        Add up over this death time, for all subjects
					*/
					ndeath =0;   /* total number of deaths at this time point */
					deathwt =0;  /* sum(wt) for the deaths */
					d_denom =0;  /*contribution to denominator for the deaths*/
					for (k=person; k< coxPh->StrataVec()[istrat]; k++)
					{
						p2 = coxPh->EndTimeIndices()[k];
						if (pData->y_ptr(1)[p2]  < dtime) break;  /* only tied times */

						nrisk++;
						denom += pData->weight_ptr()[p2] * exp(eta[p2] + pData->offset_ptr(false)[p2] - center);
						esum += eta[p2];
						if (coxPh->StatusVec()[p2] ==1)
						{
							ndeath ++;
							deathwt += pData->weight_ptr()[p2];
							d_denom += pData->weight_ptr()[p2] * exp(eta[p2] + pData->offset_ptr(false)[p2] - center);
							loglik  += pData->weight_ptr()[p2] *(eta[p2] + pData->offset_ptr(false)[p2] - center);
						}
					}
					ksave = k;

					/* compute the increment in hazard
					** hazard = usual increment
					** e_hazard = efron increment, for tied deaths only
					*/
					if (coxPh->TieApproxMethod()==0 || ndeath==1)
					{ /* Breslow */
						loglik -= deathwt*log(denom);
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
							temp = (double)k /ndeath;    /* don't do integer division*/
							loglik -= deathwt *log(denom - temp*d_denom);
							hazard += deathwt/(denom - temp*d_denom);
							e_hazard += (1-temp) *deathwt/(denom - temp*d_denom);
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
						p2 = coxPh->EndTimeIndices()[person];
						if (coxPh->StatusVec()[p2] ==1) resid[p2] = 1 + temp*exp(eta[p2] + pData->offset_ptr(false)[p2] - center);
						else resid[p2] = cumhaz * exp(eta[p2] + pData->offset_ptr(false)[p2] - center);
					}
					cumhaz += hazard;

					/* see if we need to shift the centering (very rare case) */
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
					stratastart = person;
					for (; indx1< coxPh->StrataVec()[istrat]; indx1++)
					{
						p1 = coxPh->StartTimeIndices()[indx1];
						resid[p1] -= cumhaz * exp(eta[p1] + pData->offset_ptr(false)[p1] - center);
					}
					cumhaz =0;
					denom = 0;
					istrat++;
				}
			}
	    }
	    return(loglik);
	}
};
#endif //__countingCoxState_h__
