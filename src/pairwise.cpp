// Implementation file for 'pairwise' distribution
//
// Author: Stefan Schroedl (schroedl@a9.com)

#include "pairwise.h"
#include "gbm_functions.h"
#include <limits>
#include <iostream>
#include <vector>
#include <algorithm>

void CRanker::Init(unsigned int max_items_per_group) {
  // Allocate sorting buffers
  score_rank_vec_.resize(max_items_per_group);
  ptrs_to_score_rank_vec_.resize(max_items_per_group);
}

bool CRanker::SetGroupScores(const double* const kScores,
                             const unsigned int kNumItems) {
  const double eps = 1e-10;

  if (kNumItems > score_rank_vec_.size()) {
    // Allocate additional space
    // (We should never get here if CPairwise::Initialize has been called
    // before, as expected)
    Init(kNumItems);
  }
  this->num_items_ = kNumItems;

  // Copy scores to buffer, and
  // initialize pointer array to score entries

  for (unsigned int i = 0; i < kNumItems; i++) {
    // Add small random number to break possible ties
    score_rank_vec_[i].first = kScores[i] + eps * (unif_rand() - 0.5);

    ptrs_to_score_rank_vec_[i] = &(score_rank_vec_[i]);
  }

  return true;
}

// Auxiliary struct to compare pair pointers
// decreasing order based on the first component (score)
struct CDoubleUintPairPtrComparison {
  bool operator()(const CRanker::CDoubleUintPair* kLhs,
                  const CRanker::CDoubleUintPair* kRhs) {
    return (kLhs->first > kRhs->first);
  }
};

bool CRanker::Rank() {
  // Sort the pointer array, based on decreasing score

  CDoubleUintPairPtrComparison comp;

  sort(ptrs_to_score_rank_vec_.begin(),
       ptrs_to_score_rank_vec_.begin() + num_items_, comp);

  bool is_changed = false;

  // Create inverted rank lookup

  for (unsigned int i = 0; i < num_items_; i++) {
    // Note: ranks are 1-based
    const unsigned int kNewRank = i + 1;
    if (!is_changed) {
      is_changed = (kNewRank != ptrs_to_score_rank_vec_[i]->second);
    }
    // Store the rank with the corresponding score in the vecdipScoreRank array
    ptrs_to_score_rank_vec_[i]->second = kNewRank;
  }

  return is_changed;
}

void CConc::Init(unsigned long max_group, unsigned long max_items_per_group,
                 unsigned int rank_cutoff) {
  CIRMeasure::Init(max_group, max_items_per_group, rank_cutoff);
  paircount_vec_.resize(max_group + 1, -1);
}

unsigned int CConc::PairCount(unsigned int group, const double* const kResponse,
                              unsigned int num_items) {
  if (group >= paircount_vec_.size()) {
    // Allocate additional space
    // (We should never get here if CPairwise::Initialize has been called
    // before, as expected)
    paircount_vec_.resize(group + 1, -1);
  }

  if (paircount_vec_[group] < 0.0) {
    // Not yet initialized
    paircount_vec_[group] = ComputePairCount(kResponse, num_items);
  }
  return paircount_vec_[group];
}

// Calculate the number of pairs with different labels, and store in
// veccPairCount
// Assumption: instances are sorted such that labels are non-increasing
int CConc::ComputePairCount(const double* const kResponse,
                            unsigned int num_items) {
  if (!any_pairs(kResponse, num_items)) {
    return 0;
  }

  double label_current = kResponse[0];
  int labelend = 0;  // End of range with higher labels
  int pairs = 0;

  for (unsigned int j = 1; j < num_items; j++) {
    if (kResponse[j] != label_current) {
      // i.e., dYj < dLabelCurrent
      labelend = j;
      label_current = kResponse[j];
    }
    // All items in 0 .. iLabelEnd - 1 are better than item j;
    // i.e, we have pairs (j,0), (j,1), ... (j, iLabelEnd - 1)
    pairs += labelend;
  }

  return pairs;
}

// Count the number of correctly ranked pairs with different labels
double CConc::Measure(const double* const kResponse, const CRanker& kRanker) {
  double label_current = kResponse[0];
  int label_end = 0;  // End of the range with higher labels
  int goodpairs = 0;

  for (unsigned int j = 1; j < kRanker.GetNumItems(); j++) {
    const double kYj = kResponse[j];

    if (kYj != label_current) {
      // i.e., dYj < dLabelCurrent
      label_end = j;
      label_current = kYj;
    }

    // All items in 0 .. iLabelEnd - 1 are better than this item

    for (int i = 0; i < label_end; i++) {
      if (kRanker.GetRank(i) < kRanker.GetRank(j)) {
        goodpairs++;
      }
    }
  }

  return goodpairs;
}

double CConc::SwapCost(int item_better, int item_worse,
                       const double* const kResponse,
                       const CRanker& kRanker) const {
  // Note: this implementation can handle arbitrary non-negative target values.
  // For binary (0/1) targets, the swap cost would reduce to the much simpler
  // expression:
  // (int)ranker.GetRank(iItemBetter) - (int)ranker.GetRank(iItemWorse)

  const unsigned int kRankBetter = kRanker.GetRank(item_better);
  const unsigned int kRankWorse = kRanker.GetRank(item_worse);

  // Which one of the two has the higher rank?

  unsigned int rank_upper, rank_lower;
  double resp_upper, resp_lower;
  int diff;

  if (kRankBetter > kRankWorse) {
    // Concordance increasing
    rank_upper = kRankWorse;
    rank_lower = kRankBetter;
    resp_upper = kResponse[item_worse];
    resp_lower = kResponse[item_better];

    diff = 1;  // The direct impact of the pair (iItemBetter, iItemWorse)
  } else {
    // Concordance decreasing
    rank_upper = kRankBetter;
    rank_lower = kRankWorse;
    resp_upper = kResponse[item_better];
    resp_lower = kResponse[item_worse];

    diff = -1;  // // The direct impact of the pair (iItemBetter, iItemWorse)
  }

  // Compute indirect impact for pairs involving items in between the two

  for (unsigned int rank = rank_upper + 1; rank < rank_lower; rank++) {
    const double kYi = kResponse[kRanker.GetItem(rank)];

    double score_diff = kYi - resp_lower;
    if (score_diff != 0) {
      diff += (score_diff < 0) ? 1 : -1;
    }

    score_diff = kYi - resp_upper;
    if (score_diff != 0) {
      diff += (score_diff < 0) ? -1 : 1;
    }
  }

  return diff;
}

void CNDCG::Init(unsigned long max_group, unsigned long max_items_per_group,
                 unsigned int rank_cutoff) {
  CIRMeasure::Init(max_group, max_items_per_group, rank_cutoff);

  // Initialize rank weights (note: ranks are 1-based)

  rankweight_vec_.resize(max_items_per_group + 1, 0.0);

  const unsigned int kMaxRank =
      std::min((unsigned int)max_items_per_group, get_cutoff_rank());

  // Precompute rank weights
  for (unsigned int i = 1; i <= kMaxRank; i++) {
    rankweight_vec_[i] = std::log((double)2) / log((double)(i + 1));
  }

  // Allocate buffer
  maxdcg_vec_.resize(max_group + 1, -1.0);
}

// Sum of target values, weighted by rank weight
double CNDCG::Measure(const double* const kResponse, const CRanker& kRanker) {
  double score = 0;

  for (unsigned int i = 0; i < kRanker.GetNumItems(); i++) {
    score += kResponse[i] * rankweight_vec_[kRanker.GetRank(i)];
  }
  return score;
}

double CNDCG::MaxMeasure(unsigned int group, const double* const kResponse,
                         unsigned int num_items) {
  if (group >= maxdcg_vec_.size()) {
    // Allocate additional space
    // (We should never get here if CPairwise::Initialize has been called
    // before, as expected)
    // JHickey: Then should we call initialize in the constructor?
    maxdcg_vec_.resize(group + 1, -1.0);
  }

  if (maxdcg_vec_[group] < 0.0) {
    // Not initialized

    if (!any_pairs(kResponse, num_items)) {
      // No training pairs exist
      maxdcg_vec_[group] = 0.0;
    } else {
      // Compute maximum possible DCG.
      // Note: By assumption, items are pre-sorted by descending score.

      double dScore = 0;
      unsigned int i = 0;

      while (i < num_items && kResponse[i] > 0) {
        // Note: Due to sorting, we can terminate early for a zero score.
        dScore += kResponse[i] * rankweight_vec_[i + 1];
        i++;
      }

      maxdcg_vec_[group] = dScore;
    }
  }

  return maxdcg_vec_[group];
}

double CNDCG::SwapCost(int item_better, int item_worse,
                       const double* const kResponse,
                       const CRanker& kRanker) const {
  const unsigned int kRanki = kRanker.GetRank(item_better);
  const unsigned int kRankj = kRanker.GetRank(item_worse);
  return (rankweight_vec_[kRanki] - rankweight_vec_[kRankj]) *
         (kResponse[item_better] - kResponse[item_worse]);
}

// Auxiliary function to find the top rank of a positive item (cRankTop), and
// the number of positive items (cPos)

inline void TopRankPos(const double* const kResponse, const CRanker& kRanker,
                       unsigned int& ranktop, unsigned int& pos) {
  const unsigned int kNumItems = kRanker.GetNumItems();

  ranktop = kNumItems + 1;  // Ranks are 1-based

  for (pos = 0; pos < kNumItems; pos++) {
    if (kResponse[pos] <= 0.0) {
      // All subsequent items are zero, because of presorting
      return;
    }
    ranktop = std::min(ranktop, kRanker.GetRank(pos));
  }
}

double CMRR::Measure(const double* const kResponse, const CRanker& kRanker) {
  unsigned int ranktop, pos;

  TopRankPos(kResponse, kRanker, ranktop, pos);

  const unsigned int kNumItems = std::min(kRanker.GetNumItems(), get_cutoff_rank());

  if (ranktop >= kNumItems + 1) {
    // No positive item found
    return 0.0;
  }
  // Ranks start at 1
  return 1.0 / ranktop;
}

double CMRR::SwapCost(int item_pos, int item_neg, const double* const kResponse,
                      const CRanker& kRanker) const {
  unsigned int ranktop, pos;

  TopRankPos(kResponse, kRanker, ranktop, pos);

  const unsigned int kNumItems = kRanker.GetNumItems();

  if (ranktop >= kNumItems + 1  // No positive item (ranks are 1-based)
      || pos >= kNumItems)      // No negative item
  {
    return 0.0;
  }

  const unsigned int kRankPos = kRanker.GetRank(item_pos);
  const unsigned int kRankNeg = kRanker.GetRank(item_neg);

  const unsigned int kCutoffRank = get_cutoff_rank();
  const double kMeasureCurrent = (ranktop > kCutoffRank) ? 0.0 : 1.0 / ranktop;
  const double kMeasureNeg = (kRankNeg > kCutoffRank) ? 0.0 : 1.0 / kRankNeg;

  // Only pairs where the negative item is above the top positive result,
  // or else where the positive item *is* the top item, can change the MRR

  return ((kRankNeg < ranktop || kRankPos == ranktop)
              ? (kMeasureNeg - kMeasureCurrent)
              : 0.0);
}

void CMAP::Init(unsigned long max_group, unsigned long max_items_per_group,
                unsigned int rank_cutoff) {
  CIRMeasure::Init(max_group, max_items_per_group, rank_cutoff);

  // Allocate rank buffer (note: ranks are 1-based)
  rankpos_vec_.resize(max_items_per_group + 1);
}

// Auxiliary function to find the sorted ranks of positive items (veccRankPos),
// and their number (cPos)
inline void SortRankPos(const double* const kResponse, 
                        const CRanker& kRanker,
                        std::vector<int>& rankpos_vec, 
                        unsigned int& pos) {
  // Store all ranks of positive items in veccRankPos
  for (pos = 0; pos < kRanker.GetNumItems(); pos++) {
    if (kResponse[pos] <= 0.0) {
      // All subsequent items are zero, because of presorting
      break;
    }
    rankpos_vec[pos] = kRanker.GetRank(pos);
  }

  sort(rankpos_vec.begin(), rankpos_vec.begin() + pos);
}

double CMAP::SwapCost(int item_pos, int item_neg, const double* const kResponse,
                      const CRanker& kRanker) const {
  unsigned int pos;

  SortRankPos(kResponse, kRanker, rankpos_vec_, pos);

  if (pos == 0) {
    return 0.0;
  }

  // Now veccRankPos[i] is the i-th highest rank of a positive item, and
  // cPos is the total number of positive items.

  const int kRankItemPos = kRanker.GetRank(item_pos);
  const int kRankItemNeg = kRanker.GetRank(item_neg);

  // Search for the position of the two items to swap
  const std::vector<int>::iterator kItItemPos = upper_bound(
      rankpos_vec_.begin(), rankpos_vec_.begin() + pos, kRankItemPos);
  const std::vector<int>::iterator kItItemNeg = upper_bound(
      rankpos_vec_.begin(), rankpos_vec_.begin() + pos, kRankItemNeg);

  // The number of positive items up to and including iItemPos
  const int kNumPosNotBelowItemPos = (int)(kItItemPos - rankpos_vec_.begin());

  // The number of positive items up to iItemNeg (Note: Cannot include iItemNeg
  // itself)
  const unsigned int kNumPosAboveItemNeg =
      (unsigned int)(kItItemNeg - rankpos_vec_.begin());

  // Range of indices of positive items between iRankItemPos and iRankItemNeg
  // (exclusively)
  int intermediatehigh, intermediatelow;

  // Current contribution of iItemPos
  double contrib_before = (double)kNumPosNotBelowItemPos / kRankItemPos;

  double sign, contrib_after;

  if (kRankItemNeg > kRankItemPos) {
    // MAP is decreasing
    sign = -1.0;

    // The first positive item after iRankItemPos
    intermediatelow = kNumPosNotBelowItemPos;

    // The last positive item before iRankItemNeg
    intermediatehigh = kNumPosAboveItemNeg - 1;

    // Note: iItemPos already counted in cNumPosAboveItemNeg
    contrib_after = (double)kNumPosAboveItemNeg / kRankItemNeg;
  } else {
    // MAP is increasing
    sign = 1.0;

    // The first positive result after iRankItemNeg
    intermediatelow = kNumPosAboveItemNeg;

    // The first positive result after iRankItemPos, minus iItemPos itself
    intermediatehigh = kNumPosNotBelowItemPos - 2;

    // Note: iItemPos not yet counted in cNumPosAboveItemNeg
    contrib_after = (double)(kNumPosAboveItemNeg + 1) / kRankItemNeg;
  }

  // The direct effect of switching iItemPos
  double diff = contrib_after - contrib_before;

  // The indirect effect for all items in between the two items
  for (int j = intermediatelow; j <= intermediatehigh; j++) {
    diff += sign / rankpos_vec_[j];
  }

  return diff / pos;
}

double CMAP::Measure(const double* const kResponse, const CRanker& kRanker) {
  unsigned int kPos;

  SortRankPos(kResponse, kRanker, rankpos_vec_, kPos);

  if (kPos == 0) {
    return 0.0;
  }

  // Now veccRankPos[i] is the i-th highest rank of a positive item

  double prec = 0.0;
  for (unsigned int j = 0; j < kPos; j++) {
    prec += double(j + 1) / rankpos_vec_[j];
  }

  return prec / kPos;
}

CPairwise::CPairwise(Rcpp::NumericVector misc_vec,
		     const char* kIrMeasure,
                     int num_training_rows) : misc_vec_(misc_vec) {
  // Set up adGroup - this is not required
  kGroups_ = misc_vec_.begin();

  // Set up the number of groups - this used externally
  set_num_groups(gbm_functions::NumGroups(kGroups_, num_training_rows));

  // Construct the IR Measure
  if (!strcmp(kIrMeasure, "conc")) {
    pirm_.reset(new CConc());
  } else if (!strcmp(kIrMeasure, "map")) {
    pirm_.reset(new CMAP());
  } else if (!strcmp(kIrMeasure, "mrr")) {
    pirm_.reset(new CMRR());
  } else {
    if (strcmp(kIrMeasure, "ndcg")) {
      Rprintf(
          "Unknown IR measure '%s' in initialization, using 'ndcg' instead\n",
          kIrMeasure);
    }
    pirm_.reset(new CNDCG());
  }
}

CDistribution* CPairwise::Create(DataDistParams& distparams) {
  // Create pointers to pairwise
  Rcpp::NumericVector misc_vec(distparams.misc[0]);
  
  std::size_t offset_tomeasure = distparams.family.find("_");
  if (offset_tomeasure == std::string::npos) {
    throw gbm_exception::Failure(
        "Unable to locate IR metric required for pairwise");
  }
  const char* kIrMeasure = distparams.family.c_str() + offset_tomeasure + 1;

  if (!gbm_functions::has_value(misc_vec)) {
    throw gbm_exception::Failure("Pairwise requires misc to initialize");
  } 
  return new CPairwise(misc_vec, kIrMeasure, distparams.num_trainrows);
}

CPairwise::~CPairwise() {}

// Auxiliary function for addition of optional offset parameter
inline const double* OffsetVector(const double* const kCovariates,
                                  const double* const kOffset,
                                  unsigned int start, unsigned int end,
                                  std::vector<double>& buffer_vec) {
  if (!kOffset) {
    // Optional second argument is not set, just return first one
    return kCovariates + start;
  } else {
    for (unsigned int i = start, iOut = 0; i < end; i++, iOut++) {
      buffer_vec[iOut] = kCovariates[i] + kOffset[i];
    }
    return &buffer_vec[0];
  }
}

void CPairwise::ComputeWorkingResponse(const CDataset& kData, const Bag& kBag,
                                       const double* kFuncEstimate,
                                       std::vector<double>& residuals) {
  if (kData.get_trainsize() <= 0) return;

  // Iterate through all groups, compute gradients

  unsigned int item_start = 0;
  unsigned int item_end = 0;

  while (item_start < kData.get_trainsize()) {
    residuals[item_end] = 0;
    hessian_[item_end] = 0;

    const double dGroup = kGroups_[item_start];

    // Find end of current group, initialize working response
    for (item_end = item_start + 1;
         item_end < kData.get_trainsize() && kGroups_[item_end] == dGroup;
         item_end++) {
      // Clear gradients from last iteration
      residuals[item_end] = 0;
      hessian_[item_end] = 0;
    }

    if (kBag.get_element(item_start)) {
      // Group is part of the training set

      const int kNumItems = item_end - item_start;

      // If offset given, add up current scores
      const double* kFuncPlusOffset =
          OffsetVector(kFuncEstimate, kData.offset_ptr(), item_start, item_end,
                       func_est_plus_offset_);

      // Accumulate gradients
      // TODO: Implement better way to ensure casting robust to overflow
      int int_group = 0;
      if (fabs(dGroup) > nextafter(INT_MAX, 0) || std::isnan(dGroup)) {
        int_group = copysign(INT_MAX, dGroup);
      } else {
        int_group = (int)dGroup;
      }

      ComputeLambdas(int_group, kNumItems, kData.y_ptr() + item_start,
                     kFuncPlusOffset, kData.weight_ptr() + item_start,
                     &residuals[item_start], &hessian_[item_start]);
    }

    // Next group
    item_start = item_end;
  }
}

// Referring to MSR-TR-2010-82-2, section 7 (see also the vignette):
//
// Let P be the set of pairs (i,j) where Y(i)>Y(j) (i is better than j).
// The approximation to the IR measure is the utility function C (to be
// maximized)
//   C
//   = \Sum_{(i,j) in P} |Delta Z_ij| C(s_i - s_j)
//   = \Sum_{(i,j) in P} |Delta Z_ij| / (1 + std::exp(-(s_i - s_j))),
// where |Delta Z_ij| is the cost of swapping (only) i and j in the current
// ranking,
// and s_i, s_j are the prediction scores (sum of the tree predictions) for
// items
// i and j.
//
// For (i,j) in P, define
//   lambda_ij
//   = dC(s_i-s_j) / ds_i
//   = - |Delta Z_ij| / (1 + std::exp(s_i - s_j))
//   = - |Delta Z_ij| * rho_ij,
// with
//   rho_ij = - lambda_ij / |Delta Z_ij| = 1 / (1 + std::exp(s_i - s_j))
//
// So the gradient of C with respect to s_i is
//   dC / ds_i
//   =(def) lambda_i
//   = \Sum_{j|(i,j) in P} lambda_ij - \Sum_{j|(j,i) in P} lambda_ji
//   = - \Sum_{j|(i,j) in P} |Delta Z_ij| * rho_ij
//     + \Sum_{j|(j,i) in P} |Delta Z_ji| * rho_ji;
// it is stored in adZ[i].
//
// The second derivative is
//   d^2C / ds_i^2
//   =(def) gamma_i
//   =   \Sum_{j|(i,j) in P} |Delta Z_ij| * rho_ij * (1-rho_ij)
//     - \Sum_{j|(j,i) in P} |Delta Z_ji| * rho_ji * (1-rho_ji);
// it is stored in vecdHessian[i].
//
// The Newton step for a particular leaf node is (a fraction of)
// g'/g'', where g' (resp. g'') is the sum of dC/ds_i = lambda_i
// (resp. d^2C/d^2s_i = gamma_i) over all instances falling into this leaf. This
// summation is calculated later in CPairwise::FitBestConstant().

void CPairwise::ComputeLambdas(int group, unsigned int num_items,
                               const double* const kResponse,
                               const double* const kFuncEstimate,
                               const double* const kWeights, double* residuals,
                               double* deriv) {
  // Assumption: Weights are constant within group
  if (kWeights[0] <= 0) {
    return;
  }

  // Normalize for maximum achievable group score
  const double kMaxScore = pirm_->MaxMeasure(group, kResponse, num_items);

  if (kMaxScore <= 0.0) {
    // No pairs
    return;
  }

  // Rank items by current score
  ranker_.SetGroupScores(kFuncEstimate, num_items);
  ranker_.Rank();

  double label_current = kResponse[0];

  // First index of instance that has dLabelCurrent
  // (i.e., each smaller index corresponds to better item)
  unsigned int label_current_start = 0;

  // Number of pairs with unequal labels
  unsigned int pairs = 0;

  for (unsigned int j = 1; j < num_items; j++) {
    const double kYj = kResponse[j];

    if (kYj != label_current) {
      label_current_start = j;
      label_current = kYj;
    }

    for (unsigned int i = 0; i < label_current_start; i++) {
      // Instance i is better than j

      const double kSwapCost = fabs(pirm_->SwapCost(i, j, kResponse, ranker_));

      if (!std::isfinite(kSwapCost)) {
        throw gbm_exception::Failure("infinite swap cost");
      }

      if (kSwapCost > 0.0) {
        pairs++;
        const double kRhoij =
            1.0 / (1.0 + std::exp(kFuncEstimate[i] - kFuncEstimate[j]));
        if (!std::isfinite(kRhoij)) {
          throw gbm_exception::Failure("unanticipated infinity");
        };

        const double kLambdaij = kSwapCost * kRhoij;
        residuals[i] += kLambdaij;
        residuals[j] -= kLambdaij;
        const double kDerivij = kLambdaij * (1.0 - kRhoij);
        if (kDerivij < 0) {
          throw gbm_exception::Failure("negative derivative!");
        }
        deriv[i] += kDerivij;
        deriv[j] += kDerivij;
      }
    }
  }

  if (pairs > 0) {
    // Normalize for number of training pairs
    const double kQNorm = 1.0 / (kMaxScore * pairs);

    for (unsigned int j = 0; j < num_items; j++) {
      residuals[j] *= kQNorm;
      deriv[j] *= kQNorm;
    }
  }
}

void CPairwise::Initialize(const CDataset& kData) {
  if (kData.nrow() <= 0) return;

  // Allocate memory for derivative buffer
  hessian_.resize(kData.nrow());

  // Count the groups and number of items per group
  unsigned int max_items_per_group = 0;
  double max_group = 0;

  unsigned int item_start = 0;
  unsigned int item_end = 0;

  while (item_start < kData.nrow()) {
    const double kGroup = kGroups_[item_start];

    // Find end of current group
    for (item_end = item_start + 1;
         item_end < kData.nrow() && kGroups_[item_end] == kGroup; item_end++)
      ;

    const unsigned int kNumItems = item_end - item_start;
    if (kNumItems > max_items_per_group) {
      max_items_per_group = kNumItems;
    }
    if (kGroup > max_group) {
      max_group = kGroup;
    }

    // Next group
    item_start = item_end;
  }

  // Allocate buffer for offset addition
  func_est_plus_offset_.resize(max_items_per_group);

  // Allocate ranker memory
  ranker_.Init(max_items_per_group);

  // Allocate IR measure memory

  // The last element of adGroup specifies the cutoff
  // (zero means no cutoff)
  unsigned int rank_cutoff = max_items_per_group;
  if (kGroups_[kData.nrow()] > 0) {
    rank_cutoff = (unsigned int)kGroups_[kData.nrow()];
  }

  // TODO: Make More robust against overflow
  unsigned long max_group_unsigned_long = 0;
  if (fabs(max_group) > nextafter(ULONG_MAX, 0) || std::isnan(max_group)) {
    max_group_unsigned_long = copysign(ULONG_MAX, max_group);
  } else {
    max_group_unsigned_long = (unsigned long)max_group;
  }
  pirm_->Init(max_group_unsigned_long, max_items_per_group, rank_cutoff);
}

double CPairwise::InitF(const CDataset& kData) { return 0.0; }

double CPairwise::Deviance(const CDataset& kData, const Bag& kBag,
                           const double* kfuncEstimate) {
  // Shift adGroup to validation set if necessary
  long num_rows_in_set = kData.get_size_of_set();
  if (num_rows_in_set <= 0) {
    return 0;
  }

  double loss = 0.0;
  double weight = 0.0;

  unsigned int item_start = 0;
  unsigned int item_end = item_start;
  const unsigned int kEnd = num_rows_in_set;

  while (item_start < kEnd) {
    const double kGroup = kGroups_[item_start];
    const double kWi = kData.weight_ptr()[item_start];

    // Find end of current group
    for (item_end = item_start + 1;
         item_end < kEnd && kGroups_[item_end] == kGroup; item_end++)
      ;

    const int cNumItems = item_end - item_start;
    // TODO: Implement better way to ensure casting robust to overflow
    int int_group = 0;
    if (fabs(kGroup) > nextafter(INT_MAX, 0) || std::isnan(kGroup)) {
      int_group = copysign(INT_MAX, kGroup);
    } else {
      int_group = (int)kGroup;
    }
    const double kMaxScore =
        pirm_->MaxMeasure(int_group, kData.y_ptr() + item_start, cNumItems);

    if (kMaxScore > 0.0) {
      // Rank items by current score

      // If offset given, add up current scores
      const double* kFuncPlusOffset =
          OffsetVector(kfuncEstimate, kData.offset_ptr(), item_start, item_end,
                       func_est_plus_offset_);

      ranker_.SetGroupScores(kFuncPlusOffset, cNumItems);
      ranker_.Rank();

      loss +=
          kWi * pirm_->Measure(kData.y_ptr() + item_start, ranker_) / kMaxScore;
      weight += kWi;
    }
    // Next group
    item_start = item_end;
  }

  // Loss = 1 - utility
  return 1.0 - loss / weight;
}

void CPairwise::FitBestConstant(const CDataset& kData, const Bag& kBag,
                                const double* kFuncEstimate,
                                unsigned long num_terminalnodes,
                                std::vector<double>& residuals,
                                CCARTTree& tree) {
  // Assumption: ComputeWorkingResponse() has been executed before with
  // the same arguments

  // Allocate space for numerators and denominators, and set to zero
  fit_numerator_.reserve(num_terminalnodes);
  fit_denominator_.reserve(num_terminalnodes);
  for (unsigned int i = 0; i < num_terminalnodes; i++) {
    fit_numerator_[i] = 0.0;
    fit_denominator_[i] = 0.0;
  }

  for (unsigned int obs_num = 0; obs_num < kData.get_trainsize(); obs_num++) {
    if (kBag.get_element(obs_num)) {
      fit_numerator_[tree.get_node_assignments()[obs_num]] +=
          kData.weight_ptr()[obs_num] * residuals[obs_num];
      fit_denominator_[tree.get_node_assignments()[obs_num]] +=
          kData.weight_ptr()[obs_num] * hessian_[obs_num];
    }
  }

  for (unsigned int node_num = 0; node_num < num_terminalnodes; node_num++) {
    if (tree.has_node(node_num)) {
      if (fit_denominator_[node_num] <= 0.0) {
        tree.get_terminal_nodes()[node_num]->set_prediction(0.0);
      } else {
        tree.get_terminal_nodes()[node_num]->set_prediction(
            fit_numerator_[node_num] / fit_denominator_[node_num]);
      }
    }
  }
}

double CPairwise::BagImprovement(const CDataset& kData, const Bag& kBag,
                                 const double* kFuncEstimate,
                                 const double kShrinkage,
                                 const std::vector<double>& kDeltaEstimate) {
  if (kData.get_trainsize() <= 0) {
    return 0;
  }

  double loss = 0.0;
  double weight = 0.0;

  unsigned int item_start = 0;
  unsigned int item_end = 0;

  while (item_start < kData.get_trainsize()) {
    const double kGroup = kGroups_[item_start];

    // Find end of current group
    for (item_end = item_start + 1;
         item_end < kData.get_trainsize() && kGroups_[item_end] == kGroup;
         item_end++)
      ;

    if (!kBag.get_element(item_start)) {
      // Group was held out of training set

      const unsigned int kNumItems = item_end - item_start;
      // TODO: Implement better way to ensure casting robust to overflow
      int int_group = 0;
      if (fabs(kGroup) > nextafter(INT_MAX, 0) || std::isnan(kGroup)) {
        int_group = copysign(INT_MAX, kGroup);
      } else {
        int_group = (int)kGroup;
      }
      const double dMaxScore =
          pirm_->MaxMeasure(int_group, kData.y_ptr() + item_start, kNumItems);

      if (dMaxScore > 0.0) {
        // If offset given, add up current scores
        const double* kFuncPlusOffset =
            OffsetVector(kFuncEstimate, kData.offset_ptr(), item_start,
                         item_end, func_est_plus_offset_);

        // Compute score according to old score, adF
        ranker_.SetGroupScores(kFuncPlusOffset, kNumItems);
        ranker_.Rank();
        const double kOldScore =
            pirm_->Measure(kData.y_ptr() + item_start, ranker_);

        // Compute score according to new score: adF' =  adF + dStepSize *
        // adFadj
        for (unsigned int i = 0; i < kNumItems; i++) {
          ranker_.AddToScore(i, kDeltaEstimate[i + item_start] * kShrinkage);
        }

        const double kWi = kData.weight_ptr()[item_start];

        if (ranker_.Rank()) {
          // Ranking changed
          const double kNewScore =
              pirm_->Measure(kData.y_ptr() + item_start, ranker_);
          loss += kWi * (kNewScore - kOldScore) / dMaxScore;
        }
        weight += kWi;
      }
    }

    // Next group
    item_start = item_end;
  }

  return loss / weight;
}

void CPairwise::BagData(const CDataset& kData, Bag& bag) {
  double last_group = -1;
  bool is_chosen = false;
  unsigned int bagged = 0;
  unsigned int bagged_groups = 0;
  unsigned int seen_groups = 0;
  unsigned int total_groupsinbag =
      (unsigned long)(bag.get_bagfraction() * get_num_groups());

  if (total_groupsinbag <= 0) {
    total_groupsinbag = 1;
  }

  for (unsigned long i = 0; i < kData.get_trainsize(); i++) {
    const double kGroup = kGroups_[i];

    if (kGroup != last_group) {
      if (bagged_groups >= total_groupsinbag) {
        break;
      }

      // Group changed, make a new decision
      is_chosen = (unif_rand() * (get_num_groups() - seen_groups) <
                   total_groupsinbag - bagged_groups);
      if (is_chosen) {
        bagged_groups++;
      }
      last_group = kGroup;
      seen_groups++;
    }

    if (is_chosen) {
      bag.set_element(i);
      bagged++;
    }
  }
}
