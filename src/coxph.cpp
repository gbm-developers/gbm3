//-----------------------------------
//
// File: coxph.cpp
//
// Description: Cox partial hazards model.
//
//-----------------------------------

//-----------------------------------
// Definitions
//-----------------------------------
#define frac .00000001
#define recenter 50

//-----------------------------------
// Includes
//-----------------------------------
#include "coxph.h"
#include <Rinternals.h>
#include <math.h>

//----------------------------------------
// Function Members - Private
//----------------------------------------
CCoxPH::CCoxPH(double* stats, int* sortedEnd, int* sortedSt, int* strats, bool tiedTimes):
areTiedTimes(tiedTimes), sortedEndTimes(sortedEnd), sortedStartTimes(sortedSt), strata(strats)
{
	status = stats;
	isUpdatedCoxPh = true;

}


//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CCoxPH::Create(DataDistParams& distParams)
{

	// Initialize variables to pass to constructor
	double* stat = 0;
	int* sortedSt = NULL;
	int* sortedEnd = NULL;
	bool tiedTimes = false;

	// Check that sorted arrays can be initialized
	Rcpp::NumericMatrix checkMatrix(distParams.sorted);
	if(!GBM_FUNC::has_value(checkMatrix(Rcpp::_, 0)))
	{
		throw GBM::failure("CoxPh - sort arrays have no values");
	}

	// Check if tied times or not
	Rcpp::IntegerMatrix sortMatrix(distParams.sorted);
	if(distParams.respY.ncol() > 2)
	{
		tiedTimes=true;
		stat = distParams.respY(Rcpp::_, 2).begin();
		sortedEnd = sortMatrix(Rcpp::_, 1).begin();
		sortedSt = sortMatrix(Rcpp::_, 0).begin();
	}
	else
	{
		stat = distParams.respY(Rcpp::_, 1).begin();
		sortedEnd = sortMatrix(Rcpp::_, 0).begin();
	}
	// Set up strata
	Rcpp::IntegerVector strats(distParams.strata);

	return new CCoxPH(stat, sortedEnd, sortedSt, strats.begin(), tiedTimes);
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


    // Implement first set of equations
    if(areTiedTimes)
    {
    }
    else
    {
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
				if(status[i]==1.0)
				{
					dTot += pData->weight_ptr()[i]/vecdRiskTot[i];
				}
				dF = adF[i] +  pData->offset_ptr(false)[i];
				adZ[i] = status[i] - std::exp(dF)*dTot;
			}
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

    // Set size and move to validation set if necessary
    long cLength = pData->get_trainSize();
	if(isValidationSet)
	{
		pData->shift_to_validation();
		status = shift_ptr(status, pData->get_trainSize());
		sortedEndTimes = shift_ptr(sortedEndTimes, pData->get_trainSize());
		sortedStartTimes = shift_ptr(sortedStartTimes, pData->get_trainSize());
		strata = shift_ptr(strata, pData->get_trainSize());
		cLength = pData->GetValidSize();
	}

	// Check if using Updated CoxPH or not
	if(isUpdatedCoxPh)
	{
		Rcpp::NumericVector martingaleResid(cLength, 0.0);
		// Check if there are tied times
		if(areTiedTimes)
		{

			 // Shift back for future calculations if required
			if(isValidationSet)
			{
				pData->shift_to_train();
				status = shift_ptr(status, -(pData->get_trainSize()));
				sortedEndTimes = shift_ptr(sortedEndTimes, -(pData->get_trainSize()));
				sortedStartTimes = shift_ptr(sortedStartTimes, -(pData->get_trainSize()));
				strata = shift_ptr(strata, -(pData->get_trainSize()));
			}

			return LogLikelihoodTiedTimes(cLength, pData->y_ptr(), pData->y_ptr(1),
					 	 	 	 	 	  status, pData->weight_ptr(), adF,
					 	 	 	 	 	  strata, sortedStartTimes, sortedEndTimes, martingaleResid.begin());
		}
		else
		{

			double loglik;
			loglik=  LogLikelihood(cLength, pData->y_ptr(), status,
								   pData->weight_ptr(), adF, strata,
								   sortedEndTimes, martingaleResid.begin());

			 // Shift back for future calculations if required
			if(isValidationSet)
			{
				pData->shift_to_train();
				status = shift_ptr(status, -(pData->get_trainSize()));
				sortedEndTimes = shift_ptr(sortedEndTimes, -(pData->get_trainSize()));
				sortedStartTimes = shift_ptr(sortedStartTimes, -(pData->get_trainSize()));
				strata = shift_ptr(strata, -(pData->get_trainSize()));
			}

			return loglik;
		}
	}

	// Original CoxPH implementation
    dTotalAtRisk = 0.0; 
    for(i=0; i!=cLength; i++)
    {
        dF = adF[i] +  pData->offset_ptr(false)[i];
        dTotalAtRisk += pData->weight_ptr()[i]*std::exp(dF);
        if(status[i]==1.0)
        {
            dL += pData->weight_ptr()[i]*(dF - std::log(dTotalAtRisk));
            dW += pData->weight_ptr()[i];
        }
    }

    // Shift back for future calculations if required
    if(isValidationSet)
    {
    	pData->shift_to_train();
    	status = shift_ptr(status, -(pData->get_trainSize()));
		sortedEndTimes = shift_ptr(sortedEndTimes, -(pData->get_trainSize()));
		sortedStartTimes = shift_ptr(sortedStartTimes, -(pData->get_trainSize()));
		strata = shift_ptr(strata, -(pData->get_trainSize()));
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

            if(status[i]==1.0)
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
            if(status[i]==1.0)
            {
                dReturnValue +=
                    data.weight_ptr()[i]*(shrinkage*adFadj[i] - std::log(dNum) + log(dDen));
                dW += data.weight_ptr()[i];
            }
        }
    }

    return dReturnValue/dW;
}

double CCoxPH::LogLikelihood(const int n, const double* time2, const double* status,
							const double* weight, const double* eta, const int* strata,
							const int* sort2,     double* resid, int method)
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
	center = eta[sort2[0]];

	for (person=0; person<n;)
	{
		p2 = sort2[person];
		if (status[p2] ==0)
		{
			/* add the subject to the risk set */
			resid[p2] = exp(eta[p2] - center) * cumhaz;
			nrisk++;
			denom  += weight[p2]* exp(eta[p2] - center);
			esum += eta[p2];
			person++;
		}
		else
		{
			dtime = time2[p2];  /* found a new, unique death time */
			/*
			**        Add up over this death time, for all subjects
			*/
			ndeath =0;
			deathwt =0;
			d_denom =0;  /*contribution to denominator by death at dtime */
			for (k=person; k<strata[istrat]; k++)
			{
				 p2 = sort2[k];
				 if (time2[p2]  < dtime) break;  /* only tied times */
				 nrisk++;
				 denom += weight[p2] * exp(eta[p2] - center);
				 esum += eta[p2];
			 if (status[p2] ==1)
			 {
				 ndeath ++;
				 deathwt += weight[p2];
				 d_denom += weight[p2] * exp(eta[p2] - center);
				 loglik  += weight[p2]*(eta[p2] - center);
			  }
			}
			ksave = k;

			/* compute the increment in hazard
			** hazard = usual increment
			** e_hazard = efron increment, for tied deaths only
			*/
			if (method==0 || ndeath==1)
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
					temp = (double)k /ndeath;    /* don't do integer division*/
					loglik -= deathwt * log(denom - temp*d_denom);
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
				 p2 = sort2[person];
				 if (status[p2] ==1) resid[p2] = 1 + temp*exp(eta[p2] - center);
				 else resid[p2] = cumhaz * exp(eta[p2] - center);
			}
			cumhaz += hazard;

			/* see if we need to shift the centering (a rare case) */
			if (fabs(esum/nrisk - center) > recenter)
			{
			 temp = esum/nrisk - center;
			 center += temp;
			 denom /=  exp(temp);
			}
		} // END OF ELSE PART
		/* clean up at the end of a strata */
		if (person == strata[istrat])
		{
			for (; indx1<strata[istrat]; indx1++)
			{
				 p2 = sort2[indx1];
				 resid[p2] -= cumhaz * exp(eta[p2] - center);
			}
			cumhaz =0;
			denom = 0;
			istrat++;
		}

	} // END OF LOOP OVER PERSON
	return(loglik);
}

double CCoxPH::LogLikelihoodTiedTimes(const int n, const double *time1, const double *time2,
									  const double* status, const double* weight, const double* eta,
									  const int* strata, const int* sort1, const int* sort2, double* resid, int method)
{
	int i, j, k, ksave;
	int person, p2, indx1, p1;
	int istrat;
	double cumhaz, hazard;
	int nrisk, ndeath;
	double deathwt, denom, temp, center;
	double esum, dtime, e_hazard;
	double loglik, d_denom;
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
	**     status - (cumhaz at end of their interval - cumhaz at start)*score
	*/

	istrat=0;
	indx1 =0;
	denom =0;  /* S in the math explanation */
	cumhaz =0;
	nrisk =0;   /* number at risk */
	esum =0;  /*cumulative eta, used for rescaling */
	center = eta[sort2[0]];

	for (person=0; person<n;)
	{
		p2 = sort2[person];

		if (status[p2] ==0)
		{
			 /* add the subject to the risk set */
			 resid[p2] = exp(eta[p2] - center) * cumhaz;
			 nrisk++;
			 denom  += weight[p2]* exp(eta[p2] - center);
			 esum += eta[p2];
			 person++;
		}
		else
		{
			dtime = time2[p2];  /* found a new, unique death time */
			 /*
			 ** Remove those subjects whose start time is to the right
			 **  from the risk set, and finish computation of their residual
			 */
			temp = denom;
			for (;  indx1 <person; indx1++)
			{
				 p1 = sort1[indx1];
				 if (time1[p1] < dtime) break; /* still in the risk set */
				 nrisk--;
				 resid[p1] -= cumhaz* exp(eta[p1] - center);
				 denom  -= weight[p1] * exp(eta[p1] - center);
				 esum -= eta[p1];
			}
			if (nrisk==0)
			{  /*everyone was removed! */
				denom =0;  /* remove any accumulated round off error */
				esum =0;
			}
			else if (denom/temp < frac)
			{
				 denom =0;
				 esum  =0;
				 for (; j>indx1 && j< person; j++)
				 {
					 esum += eta[sort2[j]];
				 }

				 center = esum/nrisk;
				 for (; j>indx1 && j< person; j++)
				 {
					 p2 = sort2[j];
					 denom += weight[p2] * exp(eta[p2] - center);
				 }
			}

			 /*
			 **        Add up over this death time, for all subjects
			 */
			 ndeath =0;   /* total number of deaths at this time point */
			 deathwt =0;  /* sum(wt) for the deaths */
			 d_denom =0;  /*contribution to denominator for the deaths*/
			 for (k=person; k<strata[istrat]; k++)
			 {
				 p2 = sort2[k];
				 if (time2[p2]  < dtime) break;  /* only tied times */
				 nrisk++;
				 denom += weight[p2] * exp(eta[p2] - center);
				 esum += eta[p2];
				 if (status[p2] ==1)
				 {
					 ndeath ++;
					 deathwt += weight[p2];
					 d_denom += weight[p2] * exp(eta[p2] - center);
					 loglik  += weight[p2] *(eta[p2] - center);
				 }
			 }
			 ksave = k;
			 /* compute the increment in hazard
			 ** hazard = usual increment
			 ** e_hazard = efron increment, for tied deaths only
			 */
			 if (method==0 || ndeath==1)
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
				 p2 = sort2[person];
				 if (status[p2] ==1) resid[p2] = 1 + temp*exp(eta[p2] - center);
				 else resid[p2] = cumhaz * exp(eta[p2] - center);
			 }
			 cumhaz += hazard;

			 /* see if we need to shift the centering (very rare case) */
			 if (fabs(esum/nrisk - center) > recenter)
			 {
				 temp = esum/nrisk - center;
				 center += temp;
				 denom /=  exp(temp);
			 }
		} // END OF ELSE OVER PERSON

		/* clean up at the end of a strata */
		if (person == strata[istrat])
		{
			 for (; indx1<strata[istrat]; indx1++)
			 {
				 p1 = sort1[indx1];
				 resid[p1] -= cumhaz * exp(eta[p1] - center);
			 }
		 cumhaz =0;
		 denom = 0;
		 istrat++;
		}
	} // END OF FOR LOOP OVER PERSON

	return(loglik);
}



