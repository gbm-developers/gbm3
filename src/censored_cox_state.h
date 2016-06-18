//------------------------------------------------------------------------------
//
//  File:       censoredCoxState.h
//
//  Description: censored CoxPH methods
//
//	Author: 	James Hickey
//------------------------------------------------------------------------------

#ifndef CENSOREDCOXSTATE_H
#define CENSOREDCOXSTATE_H

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
class CensoredCoxState : public GenericCoxState {
 public:
  //----------------------
  // Public Constructors
  //----------------------
  CensoredCoxState(CCoxPH* coxPhPtr) : coxph_(coxPhPtr){};

  //---------------------
  // Public destructor
  //---------------------
  ~CensoredCoxState() {};

  //---------------------
  // Public Functions
  //---------------------
  void ComputeWorkingResponse(const CDataset& kData, const Bag& kBag,
                              const double* kFuncEstimate,
                              std::vector<double>& residuals) {
    // Initialize parameters
    std::vector<double> martingale_resid(kData.get_trainsize(), 0.0);
    LogLikelihood(kData.get_trainsize(), kData, kBag, kFuncEstimate,
                  &martingale_resid[0], false);

    // Fill up response
    for (unsigned long i = 0; i < kData.get_trainsize(); i++) {
      if (kBag.get_element(i)) {
        residuals[i] =
            kData.weight_ptr()[i] * martingale_resid[i];  // From chain rule
      }
    }
  }

  void FitBestConstant(const CDataset& kData, const Bag& kBag,
                       const double* kFuncEstimate,
                       unsigned long num_terminalnodes,
                       std::vector<double>& residuals, CCARTTree& tree) {
    // Calculate the expected number of events and actual number of events in
    // terminal nodes
    std::vector<double> martingale_resid(kData.get_trainsize(), 0.0);
    std::vector<double> expnum_events_in_nodes(num_terminalnodes,
                                               1.0 / coxph_->PriorCoeffVar());
    std::vector<double> num_events_in_nodes(num_terminalnodes,
                                            1.0 / coxph_->PriorCoeffVar());
    LogLikelihood(kData.get_trainsize(), kData, kBag, kFuncEstimate,
                  &martingale_resid[0], false);

    for (unsigned long i = 0; i < kData.get_trainsize(); i++) {
      if (kBag.get_element(i) &&
          (tree.get_terminal_nodes()[tree.get_node_assignments()[i]]
               ->get_numobs() >= tree.min_num_obs_required())) {
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
                  const Bag& kBag, const double* kFuncEstimate) {
    // Initialize Parameters
    double loglik = 0.0;
    std::vector<double> martingale_resid(kNumRowsInSet, 0.0);

    // Calculate Deviance - skip bag as pointing at validation set
    loglik = LogLikelihood(kNumRowsInSet, kData, kBag, kFuncEstimate,
                           &martingale_resid[0]);

    return -loglik;
  }

  double BagImprovement(const CDataset& kData, const Bag& kBag,
                        const double* kFuncEstimate, const double kShrinkage,
                        const std::vector<double>& kDeltaEstimate) {
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
        LogLikelihood(kData.get_trainsize(), kData, kBag, kFuncEstimate,
                      &martingale_resid_no_adj[0], false, false);
    loglike_with_adj =
        LogLikelihood(kData.get_trainsize(), kData, kBag, &eta_adj[0],
                      &martingale_resid_with_adj[0], false, false);

    return (loglike_with_adj - loglike_no_adj);
  }

 private:
  CCoxPH* coxph_;

  double LogLikelihood(const int n, const CDataset& kData, const Bag& kBag,
                       const double* eta, double* resid, bool skipbag = true,
                       bool checkinbag = true) {
    int k, ksave;
    int person, p2;
    int istrat, indx1;  // this counts up over the strata
    double cumhaz, hazard;
    int nrisk, ndeath;
    double deathwt, denom, temp, center;
    double esum, dtime, e_hazard;
    double loglik, d_denom;

    //  'person' walks through the the data from 1 to n, p2= sort2[person].
    //     sort2[0] points to the largest stop time, sort2[1] the next, ...
    //  'dtime' is a scratch variable holding the time of current interest

    istrat = 0;
    indx1 = 0;
    denom = 0;  // S in the math explanation
    cumhaz = 0;
    nrisk = 0;  // number at risk
    esum = 0;   // cumulative eta, used for rescaling
    loglik = 0;

    // Center based on last person in the set of interest
    double new_center = 0.0;
    center = -10.0E16;

    for (person = 0; person < n; person++) {
      p2 = coxph_->EndTimeIndices()[person];
      if (skipbag || (kBag.get_element(p2) == checkinbag)) {
        new_center = eta[coxph_->EndTimeIndices()[p2]] +
                     kData.offset_ptr()[coxph_->EndTimeIndices()[p2]];
        if (new_center > center) {
          center = new_center;
        }
      }
    }

    // Loop over patients
    for (person = 0; person < n;) {
      p2 = coxph_->EndTimeIndices()[person];

      // Check if bagging is required - p2 gives the within strata order
      if (skipbag || (kBag.get_element(p2) == checkinbag)) {
        if (coxph_->StatusVec()[p2] == 0) {
          /* add the subject to the risk set */
          resid[p2] = exp(eta[p2] + kData.offset_ptr()[p2] - center) * cumhaz;
          nrisk++;
          denom += kData.weight_ptr()[p2] *
                   exp(eta[p2] + kData.offset_ptr()[p2] - center);
          esum += eta[p2] + kData.offset_ptr()[p2];
          person++;
        } else {
          dtime = kData.y_ptr()[p2];  // found a new, unique death time

          //        Add up over this death time, for all subjects
          ndeath = 0;
          deathwt = 0;
          d_denom = 0;  // contribution to denominator by death at dtime
          for (k = person; k < coxph_->StrataVec()[istrat]; k++) {
            p2 = coxph_->EndTimeIndices()[k];
            // Check in loop over stratum that person in stratum has correct bag
            // properties
            if (skipbag || (kBag.get_element(p2) == checkinbag)) {
              if (kData.y_ptr()[p2] < dtime) {
                break;  // only tied times
              }
              nrisk++;
              denom += kData.weight_ptr()[p2] *
                       exp(eta[p2] + kData.offset_ptr()[p2] - center);
              esum += eta[p2] + kData.offset_ptr()[p2];
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
          if (coxph_->TieApproxMethod() == 0 || ndeath == 1) { /* Breslow */
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
          // occur.

          temp = cumhaz + (hazard - e_hazard);
          for (; person < ksave; person++) {
            p2 = coxph_->EndTimeIndices()[person];
            // Check if person in stratum in/out of bag
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

          // see if we need to shift the centering (a rare case)
          if (fabs(esum / nrisk - center) > recenter) {
            temp = esum / nrisk - center;
            center += temp;
            denom /= exp(temp);
          }
        }

        // clean up at the end of a strata
        if (person == coxph_->StrataVec()[istrat]) {
          for (; indx1 < coxph_->StrataVec()[istrat]; indx1++) {
            p2 = coxph_->EndTimeIndices()[indx1];

            // Check bagging status
            if (skipbag || (kBag.get_element(p2) == checkinbag)) {
              resid[p2] -=
                  cumhaz * exp(eta[p2] + kData.offset_ptr()[p2] - center);
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
#endif  // CENSOREDCOXSTATE_H
