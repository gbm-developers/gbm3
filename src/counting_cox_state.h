//------------------------------------------------------------------------------
//
//  File:       countingCoxState.h
//
//  Description: counting CoxPH methods
//
//	Author: 	James Hickey
//------------------------------------------------------------------------------

#ifndef COUNTINGCOXSTATE_H
#define COUNTINGCOXSTATE_H

//------------------------------
// Includes
//------------------------------
#include "dataset.h"
#include "distribution.h"
#include "generic_cox_state.h"
#include <Rcpp.h>

//------------------------------
// Class Definition
//------------------------------
class CountingCoxState : public GenericCoxState {
 public:
  //----------------------
  // Public Constructors
  //----------------------
  CountingCoxState(CCoxPH* coxPhPtr) : coxph_(coxPhPtr){};

  //---------------------
  // Public destructor
  //---------------------
  ~CountingCoxState() { coxph_ = NULL; };

  //---------------------
  // Public Functions
  //---------------------
  void ComputeWorkingResponse(const CDataset& kData,
		  	  	  	  	  	  const Bag& kBag,
                              const double* kFuncEstimate, std::vector<double>& residuals) {
    // Initialize parameters
    std::vector<double> martingale_resid(kData.get_trainsize(), 0.0);
    LogLikelihoodTiedTimes(kData.get_trainsize(), kData, kBag, kFuncEstimate,
                           &martingale_resid[0], false);

    // Fill up response
    for (unsigned long i = 0; i < kData.get_trainsize(); i++) {
      if (kBag.get_element(i)) {
        residuals[i] =
            kData.weight_ptr()[i] * martingale_resid[i];  // From chain rule
      }
    }
  }

  void FitBestConstant(const CDataset& kData, const Bag& kBag, const double* kFuncEstimate,
                       unsigned long num_terminalnodes, std::vector<double>& residuals,
                       CCARTTree& tree) {
    // Calculate the expected number of events and actual number of events in
    // terminal nodes
    std::vector<double> martingale_resid(kData.get_trainsize(), 0.0);
    std::vector<double> expnum_events_in_nodes(num_terminalnodes,
                                               1.0 / coxph_->PriorCoeffVar());
    std::vector<double> num_events_in_nodes(num_terminalnodes,
                                            1.0 / coxph_->PriorCoeffVar());
    LogLikelihoodTiedTimes(kData.get_trainsize(), kData, kBag, kFuncEstimate,
                           &martingale_resid[0], false);

    for (unsigned long i = 0; i < kData.get_trainsize(); i++) {
      if (kBag.get_element(i) &&
          (tree.get_terminal_nodes()[tree.get_node_assignments()[i]]->get_numobs() >=
           tree.min_num_obs_required())) {
        // Cap expected number of events to be at least 0
        expnum_events_in_nodes[tree.get_node_assignments()[i]] +=
            max(0.0, coxph_->StatusVec()[i] - martingale_resid[i]);
        num_events_in_nodes[tree.get_node_assignments()[i]] +=
            coxph_->StatusVec()[i];
      }
    }

    // Update Node predictions
    for (unsigned long nodeNum = 0; nodeNum < num_terminalnodes; nodeNum++) {
      // If there are no data points in node this is 0.0
      tree.get_terminal_nodes()[nodeNum]->set_prediction(
          log(num_events_in_nodes[nodeNum] / expnum_events_in_nodes[nodeNum]));
    }
  }

  double Deviance(const long kNumRowsInSet, const CDataset& kData,
		  	  	  const Bag& kBag,
                  const double* kFuncEstimate) {
    // Initialize Parameters
    double loglik = 0.0;
    std::vector<double> martingale_resid(kNumRowsInSet, 0.0);

    // Calculate Deviance
    loglik = LogLikelihoodTiedTimes(kNumRowsInSet, kData, kBag, kFuncEstimate,
                                    &martingale_resid[0]);

    return -loglik;
  }

  double BagImprovement(const CDataset& kData, const Bag& kBag, const double* kFuncEstimate,
                        const double kShrinkage, const std::vector<double>& kDeltaEstimate) {
    // Initialize Parameters
    double loglike_no_adj = 0.0;
    double loglike_with_adj = 0.0;
    std::vector<double> martingale_resid_no_adj(kData.get_trainsize(), 0.0);
    std::vector<double> martingale_resid_with_adj(kData.get_trainsize(), 0.0);
    std::vector<double> eta_adj(kData.get_trainsize(), 0.0);

    // Fill up the adjusted and shrunk eta
    for (unsigned long i = 0; i < kData.get_trainsize(); i++) {
      if (!kBag.get_element(i)) {
        eta_adj[i] = kFuncEstimate[i] + kShrinkage * kDeltaEstimate[i];

      } else {
        eta_adj[i] = kFuncEstimate[i];
      }
    }

    // Calculate likelihoods - data not in bags
    loglike_no_adj =
        LogLikelihoodTiedTimes(kData.get_trainsize(), kData, kBag, kFuncEstimate,
                               &martingale_resid_no_adj[0], false, false);
    loglike_with_adj =
        LogLikelihoodTiedTimes(kData.get_trainsize(), kData, kBag, &eta_adj[0],
                               &martingale_resid_with_adj[0], false, false);

    return (loglike_with_adj - loglike_no_adj);
  }

 private:
  CCoxPH* coxph_;
  double LogLikelihoodTiedTimes(const int n, const CDataset& kData,
		  	  	  	  	  	  	const Bag& kBag,
                                const double* eta, double* resid,
                                bool skipbag = true, bool checkinbag = true) {
    int k, ksave;
    int person, p2, indx1, p1;
    int istrat;
    double cumhaz, hazard;
    int nrisk, ndeath;
    double deathwt, denom, temp, center;
    double esum, dtime, e_hazard;
    double loglik, d_denom;

    //  'person' walks through the the data from 1 to n, p2= sort2[person].
    //     sort2[0] points to the largest stop time, sort2[1] the next, ...
    //  'dtime' is a scratch variable holding the time of current interest
    //  'indx1' walks through the start times.  It will be smaller than
    //    'person': if person=27 that means that 27 subjects have time2 >=dtime,
    //    and are thus potential members of the risk set.  If 'indx1' =9,
    //    that means that 9 subjects have start >=time and thus are NOT part
    //    of the risk set.  (stop > start for each subject guarrantees that
    //    the 9 are a subset of the 27). p1 = sort1[indx1]
    //  Basic algorithm: move 'person' forward, adding the new subject into
    //    the risk set.  If this is a new, unique death time, take selected
    //    old obs out of the sums, add in obs tied at this time, then update
    //    the cumulative hazard. Everything resets at the end of a stratum.
    //  The sort order is from large time to small, so we encounter a subjects
    //    ending time first, then their start time.
    //  The martingale residual for a subject is
    //     status - (cumhaz at entry - cumhaz at exit)*score

    istrat = 0;
    indx1 = 0;
    denom = 0; /* S in the math explanation */
    cumhaz = 0;
    nrisk = 0; /* number at risk */
    esum = 0;  /*cumulative eta, used for rescaling */
    loglik = 0;

    // Center based on last person in the set of interest
    double newcenter = 0.0;
    center = -10.0E16;

    for (person = 0; person < n; person++) {
      p2 = coxph_->EndTimeIndices()[person];
      if (skipbag || (kBag.get_element(p2) == checkinbag)) {
        newcenter = eta[coxph_->EndTimeIndices()[p2]] +
                    kData.offset_ptr()[coxph_->EndTimeIndices()[p2]];
        if (newcenter > center) {
          center = newcenter;
        }
      }
    }

    for (person = 0; person < n;) {
      p2 = coxph_->EndTimeIndices()[person];

      // Check if bagging is required
      if (skipbag || (kBag.get_element(p2) == checkinbag)) {
        if (coxph_->StatusVec()[p2] == 0) {
          // add the subject to the risk set
          resid[p2] = exp(eta[p2] + kData.offset_ptr()[p2] - center) * cumhaz;
          nrisk++;
          denom += kData.weight_ptr()[p2] *
                   exp(eta[p2] + kData.offset_ptr()[p2] - center);
          esum += eta[p2] + kData.offset_ptr()[p2];
          person++;
        } else {
          dtime = kData.y_ptr(1)[p2];  // found a new, unique death time

          // Remove those subjects whose start time is to the right
          //  from the risk set, and finish computation of their residual

          temp = denom;
          for (; indx1 < person; indx1++) {
            p1 = coxph_->StartTimeIndices()[indx1];
            if (skipbag || (kBag.get_element(p1) == checkinbag)) {
              if (kData.y_ptr(0)[p1] < dtime) break; /* still in the risk set */
              nrisk--;
              resid[p1] -=
                  cumhaz * exp(eta[p1] + kData.offset_ptr()[p1] - center);
              denom -= kData.weight_ptr()[p1] *
                       exp(eta[p1] + kData.offset_ptr()[p1] - center);
              esum -= eta[p1] + kData.offset_ptr()[p1];
            }
          }
          if (nrisk == 0) {
            // everyone was removed!
            // This happens with manufactured start/stop
            //  data sets that have some g(time) as a covariate.
            //  Just like a strata, reset the sums.
            denom = 0;
            esum = 0;
          }

          //        Add up over this death time, for all subjects
          ndeath = 0;   // total number of deaths at this time point
          deathwt = 0;  // sum(wt) for the deaths
          d_denom = 0;  // contribution to denominator for the deaths
          for (k = person; k < coxph_->StrataVec()[istrat]; k++) {
            p2 = coxph_->EndTimeIndices()[k];
            if (skipbag || (kBag.get_element(p2) == checkinbag)) {
              if (kData.y_ptr(1)[p2] < dtime) break;  // only tied times
              nrisk++;
              denom += kData.weight_ptr()[p2] *
                       exp(eta[p2] + kData.offset_ptr()[p2] - center);
              esum += eta[p2];
              if (coxph_->StatusVec()[p2] == 1) {
                ndeath++;
                deathwt += kData.weight_ptr()[p2];
                d_denom += kData.weight_ptr()[p2] *
                           exp(eta[p2] + kData.offset_ptr()[p2] - center);
                loglik += kData.weight_ptr()[p2] *
                          (eta[p2] + kData.offset_ptr()[p2] - center);
              }
            }
          }
          ksave = k;

          // compute the increment in hazard
          // hazard = usual increment
          // e_hazard = efron increment, for tied deaths only
          if (coxph_->TieApproxMethod() == 0 || ndeath == 1) {  // Breslow
            loglik -= deathwt * log(denom);
            hazard = deathwt / denom;
            e_hazard = hazard;
          } else {  // Efron
            hazard = 0;
            e_hazard = 0;       // hazard experienced by a tied death
            deathwt /= ndeath;  // average weight of each death
            for (k = 0; k < ndeath; k++) {
              temp = (double)k / ndeath;  // don't do integer division
              loglik -= deathwt * log(denom - temp * d_denom);
              hazard += deathwt / (denom - temp * d_denom);
              e_hazard += (1 - temp) * deathwt / (denom - temp * d_denom);
            }
          }

          // Give initial value to all intervals ending at this time
          // If tied censors are sorted before deaths (which at least some
          //  callers of this routine do), then the else below will never
          //  occur.

          temp = cumhaz + (hazard - e_hazard);
          for (; person < ksave; person++) {
            p2 = coxph_->EndTimeIndices()[person];
            if (skipbag || (kBag.get_element(p2) == checkinbag)) {
              if (coxph_->StatusVec()[p2] == 1)
                resid[p2] =
                    1 + temp * exp(eta[p2] + kData.offset_ptr()[p2] - center);
              else
                resid[p2] =
                    cumhaz * exp(eta[p2] + kData.offset_ptr()[p2] - center);
            }
          }
          cumhaz += hazard;

          // see if we need to shift the centering (very rare case)
          if (fabs(esum / nrisk - center) > recenter) {
            temp = esum / nrisk - center;
            center += temp;
            denom /= exp(temp);
          }
        }

        // clean up at the end of a strata
        if (person == coxph_->StrataVec()[istrat]) {
          for (; indx1 < coxph_->StrataVec()[istrat]; indx1++) {
            p1 = coxph_->StartTimeIndices()[indx1];
            if (skipbag || (kBag.get_element(p1) == checkinbag)) {
              resid[p1] -=
                  cumhaz * exp(eta[p1] + kData.offset_ptr()[p1] - center);
            }
          }
          cumhaz = 0;
          denom = 0;
          istrat++;
        }
      } else {
        // Increment person if not in bag
        person++;
      }
    }
    return (loglik);
  }
};
#endif  // COUNTINGCOXSTATE_H
