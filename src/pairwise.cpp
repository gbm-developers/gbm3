// Implementation file for 'pairwise' distribution
//
// Author: Stefan Schroedl (schroedl@a9.com)

#include "pairwise.h"
#include "gbm_functions.h"
#include <limits>
#include <iostream>
#include <vector>
#include <algorithm>

//#define NOISY_DEBUG
#ifdef NOISY_DEBUG
#endif

void CRanker::Init(unsigned int cMaxItemsPerGroup)
{
    // Allocate sorting buffers
    score_rank_vec_.resize(cMaxItemsPerGroup);
    ptrs_to_score_rank_vec_.resize(cMaxItemsPerGroup);
}

bool CRanker::SetGroupScores(const double* const adScores, const unsigned int cNumItems)
{
    const double dEPS = 1e-10;

    if (cNumItems > score_rank_vec_.size())
    {
        // Allocate additional space
        // (We should never get here if CPairwise::Initialize has been called before, as expected)
        Init(cNumItems);
    }
    this->num_items_ = cNumItems;

    // Copy scores to buffer, and
    // initialize pointer array to score entries

    for(unsigned int i = 0; i < cNumItems; i++)
    {
        // Add small random number to break possible ties
        score_rank_vec_[i].first = adScores[i] + dEPS * (unif_rand() - 0.5);

        ptrs_to_score_rank_vec_[i] = &(score_rank_vec_[i]);
    }

    return true;
}

// Auxiliary struct to compare pair pointers
// decreasing order based on the first component (score)
struct CDoubleUintPairPtrComparison
{
    bool operator() (const CRanker::CDoubleUintPair* lhs, const CRanker::CDoubleUintPair* rhs)
    {
        return (lhs->first > rhs->first);
    }
};

bool CRanker::Rank()
{
    // Sort the pointer array, based on decreasing score

    CDoubleUintPairPtrComparison comp;

    sort(ptrs_to_score_rank_vec_.begin(), ptrs_to_score_rank_vec_.begin() + num_items_, comp);

    bool bChanged = false;

    // Create inverted rank lookup

    for(unsigned int i = 0; i < num_items_; i++)
    {
        // Note: ranks are 1-based
        const unsigned int cNewRank = i + 1;
        if (!bChanged)
        {
            bChanged = (cNewRank != ptrs_to_score_rank_vec_[i]->second);
        }
        // Store the rank with the corresponding score in the vecdipScoreRank array
        ptrs_to_score_rank_vec_[i]->second = cNewRank;
    }

    return bChanged;
}



void CConc::Init
(
    unsigned long cMaxGroup,
    unsigned long cMaxItemsPerGroup,
    unsigned int cRankCutoff
)
{
    CIRMeasure::Init(cMaxGroup, cMaxItemsPerGroup, cRankCutoff);
    paircount_vec_.resize(cMaxGroup + 1, -1);
}

unsigned int CConc::PairCount(unsigned int iGroup, const double* const adY, unsigned int cNumItems)
{
    if (iGroup >= paircount_vec_.size())
    {
        // Allocate additional space
        // (We should never get here if CPairwise::Initialize has been called before, as expected)
        paircount_vec_.resize(iGroup + 1, -1);
    }

    if (paircount_vec_[iGroup] < 0.0)
    {
        // Not yet initialized
        paircount_vec_[iGroup] = ComputePairCount(adY, cNumItems);
    }
    return paircount_vec_[iGroup];
}

// Calculate the number of pairs with different labels, and store in veccPairCount
// Assumption: instances are sorted such that labels are non-increasing
int CConc::ComputePairCount(const double* const adY, unsigned int cNumItems)
{
    if (!any_pairs(adY, cNumItems))
    {
        return 0;
    }

    double dLabelCurrent = adY[0];
    int iLabelEnd        = 0; // End of range with higher labels
    int cPairs           = 0;

    for (unsigned int j = 1; j < cNumItems; j++)
    {
        if (adY[j] != dLabelCurrent)
        {
            // i.e., dYj < dLabelCurrent
            iLabelEnd     = j;
            dLabelCurrent = adY[j];
        }
        // All items in 0 .. iLabelEnd - 1 are better than item j;
        // i.e, we have pairs (j,0), (j,1), ... (j, iLabelEnd - 1)
        cPairs += iLabelEnd;
    }

    return cPairs;
}

// Count the number of correctly ranked pairs with different labels
double CConc::Measure(const double* const adY, const CRanker& ranker)
{
    double dLabelCurrent  = adY[0];
    int iLabelEnd         = 0; // End of the range with higher labels
    int cGoodPairs        = 0;

    for (unsigned int j = 1; j < ranker.GetNumItems(); j++)
    {
        const double dYj  = adY[j];

        if (dYj != dLabelCurrent)
        {
            // i.e., dYj < dLabelCurrent
            iLabelEnd     = j;
            dLabelCurrent = dYj;
        }

        // All items in 0 .. iLabelEnd - 1 are better than this item

        for (int i = 0; i < iLabelEnd; i++)
        {
            if (ranker.GetRank(i) < ranker.GetRank(j))
            {
                cGoodPairs++;
            }
        }
    }

    return cGoodPairs;
}

double CConc::SwapCost(int iItemBetter, int iItemWorse, const double* const adY, const CRanker& ranker) const
{
    // Note: this implementation can handle arbitrary non-negative target values.
    // For binary (0/1) targets, the swap cost would reduce to the much simpler expression:
    // (int)ranker.GetRank(iItemBetter) - (int)ranker.GetRank(iItemWorse)

    const unsigned int cRankBetter = ranker.GetRank(iItemBetter);
    const unsigned int cRankWorse  = ranker.GetRank(iItemWorse);

    // Which one of the two has the higher rank?

    unsigned int cRankUpper, cRankLower;
    double dYUpper, dYLower;
    int cDiff;

    if (cRankBetter > cRankWorse)
    {
        // Concordance increasing
        cRankUpper = cRankWorse;
        cRankLower = cRankBetter;
        dYUpper    = adY[iItemWorse];
        dYLower    = adY[iItemBetter];

        cDiff = 1; // The direct impact of the pair (iItemBetter, iItemWorse)
    }
    else
    {
        // Concordance decreasing
        cRankUpper = cRankBetter;
        cRankLower = cRankWorse;
        dYUpper    = adY[iItemBetter];
        dYLower    = adY[iItemWorse];

        cDiff = -1; // // The direct impact of the pair (iItemBetter, iItemWorse)
    }

    // Compute indirect impact for pairs involving items in between the two

    for (unsigned int cRank = cRankUpper + 1; cRank < cRankLower; cRank++)
    {
        const double dYi = adY[ranker.GetItem(cRank)];

        double dScoreDiff = dYi - dYLower;
        if (dScoreDiff != 0)
        {
            cDiff += (dScoreDiff < 0) ? 1 : -1;
        }

        dScoreDiff = dYi - dYUpper;
        if (dScoreDiff != 0)
        {
            cDiff += (dScoreDiff < 0) ? -1 : 1;
        }
    }

    return cDiff;
}


void CNDCG::Init
(
    unsigned long cMaxGroup,
    unsigned long cMaxItemsPerGroup,
    unsigned int cRankCutoff
)
{
    CIRMeasure::Init(cMaxGroup, cMaxItemsPerGroup, cRankCutoff);

    // Initialize rank weights (note: ranks are 1-based)

    rankweight_vec_.resize(cMaxItemsPerGroup + 1, 0.0);

    const unsigned int cMaxRank = std::min((unsigned int)cMaxItemsPerGroup, get_cutoff_rank());

    // Precompute rank weights
    for (unsigned int i = 1; i <= cMaxRank; i++)
    {
        rankweight_vec_[i] = std::log((double)2) / log((double)(i+1));
    }

    // Allocate buffer
    maxdcg_vec_.resize(cMaxGroup + 1, -1.0);
}


// Sum of target values, weighted by rank weight
double CNDCG::Measure(const double* const adY, const CRanker& ranker)
{
    double dScore = 0;

    for (unsigned int i = 0; i < ranker.GetNumItems(); i++)
    {
        dScore += adY[i] * rankweight_vec_[ranker.GetRank(i)];
    }
    return dScore;
}

double CNDCG::MaxMeasure(unsigned int iGroup, const double* const adY, unsigned int cNumItems)
{
    if (iGroup >= maxdcg_vec_.size())
    {
        // Allocate additional space
        // (We should never get here if CPairwise::Initialize has been called before, as expected)
    	// JHickey: Then should we call initialize in the constructor?
        maxdcg_vec_.resize(iGroup + 1, -1.0);
    }

    if (maxdcg_vec_[iGroup] < 0.0)
    {
        // Not initialized

        if (!any_pairs(adY, cNumItems))
        {
            // No training pairs exist
            maxdcg_vec_[iGroup] = 0.0;
        }
        else
        {
            // Compute maximum possible DCG.
            // Note: By assumption, items are pre-sorted by descending score.

            double dScore = 0;
            unsigned int i = 0;

            while (i < cNumItems && adY[i] > 0)
            {
                // Note: Due to sorting, we can terminate early for a zero score.
                dScore += adY[i] * rankweight_vec_[i + 1];
                i++;
            }

            maxdcg_vec_[iGroup] = dScore;

#ifdef NOISY_DEBUG
            if (maxdcg_vec_[iGroup] == 0)
            {
                Rprintf("max score is 0: iGroup = %d, maxScore = %f\n", 
                        iGroup,  maxdcg_vec_[iGroup]);
                throw GBM::Failure();
            }
#endif
        }
    }

    return maxdcg_vec_[iGroup];
}


double CNDCG::SwapCost(int iItemBetter, int iItemWorse, const double* const adY, const CRanker& ranker) const
{
    const unsigned int cRanki = ranker.GetRank(iItemBetter);
    const unsigned int cRankj = ranker.GetRank(iItemWorse);
    return (rankweight_vec_[cRanki] - rankweight_vec_[cRankj]) * (adY[iItemBetter] - adY[iItemWorse]);
}


// Auxiliary function to find the top rank of a positive item (cRankTop), and the number of positive items (cPos)

inline void TopRankPos(const double* const adY, const CRanker& ranker, unsigned int& cRankTop, unsigned int& cPos)
{
    const unsigned int cNumItems = ranker.GetNumItems();

    cRankTop = cNumItems + 1; // Ranks are 1-based

    for (cPos = 0; cPos < cNumItems; cPos++)
    {
        if (adY[cPos] <= 0.0)
        {
            // All subsequent items are zero, because of presorting
            return;
        }
        cRankTop = min(cRankTop, ranker.GetRank(cPos));
    }
}

double CMRR::Measure(const double* const adY, const CRanker& ranker)
{
    unsigned int cRankTop, cPos;

    TopRankPos(adY, ranker, cRankTop, cPos);

    const unsigned int cNumItems =  min(ranker.GetNumItems(), get_cutoff_rank());

    if (cRankTop >= cNumItems + 1)
    {
        // No positive item found
        return 0.0;
    }
    // Ranks start at 1
    return 1.0 / cRankTop;
}

double CMRR::SwapCost(int iItemPos, int iItemNeg, const double* const adY, const CRanker& ranker) const
{
    unsigned int cRankTop, cPos;

    TopRankPos(adY, ranker, cRankTop, cPos);

    const unsigned int cNumItems = ranker.GetNumItems();

    if (cRankTop >= cNumItems + 1 // No positive item (ranks are 1-based)
        || cPos >= cNumItems)     // No negative item
    {
        return 0.0;
    }

    const unsigned int cRankPos    = ranker.GetRank(iItemPos);
    const unsigned int cRankNeg    = ranker.GetRank(iItemNeg);

    const unsigned int cCutoffRank = get_cutoff_rank();
    const double dMeasureCurrent   = (cRankTop > cCutoffRank) ? 0.0 : 1.0 / cRankTop;
    const double dMeasureNeg       = (cRankNeg > cCutoffRank) ? 0.0 : 1.0 / cRankNeg;

    // Only pairs where the negative item is above the top positive result,
    // or else where the positive item *is* the top item, can change the MRR

    return ((cRankNeg < cRankTop || cRankPos == cRankTop) ? (dMeasureNeg - dMeasureCurrent) : 0.0);
}

void CMAP::Init
(
    unsigned long cMaxGroup,
    unsigned long cMaxItemsPerGroup,
    unsigned int cRankCutoff
)
{
    CIRMeasure::Init(cMaxGroup, cMaxItemsPerGroup, cRankCutoff);

    // Allocate rank buffer (note: ranks are 1-based)
    rankpos_vec_.resize(cMaxItemsPerGroup + 1);
}



// Auxiliary function to find the sorted ranks of positive items (veccRankPos), and their number (cPos)
inline void SortRankPos(const double* const adY, const CRanker& ranker, vector<int>& veccRankPos, unsigned int& cPos)
{
    // Store all ranks of positive items in veccRankPos
    for (cPos = 0; cPos < ranker.GetNumItems(); cPos++)
    {
        if (adY[cPos] <= 0.0)
        {
            // All subsequent items are zero, because of presorting
            break;
        }
        veccRankPos[cPos] = ranker.GetRank(cPos);
    }

    sort(veccRankPos.begin(), veccRankPos.begin() + cPos);
}


double CMAP::SwapCost(int iItemPos, int iItemNeg, const double* const adY, const CRanker& ranker) const
{
    unsigned int cPos;

    SortRankPos(adY, ranker, rankpos_vec_, cPos);

    if (cPos == 0)
    {
        return 0.0;
    }

    // Now veccRankPos[i] is the i-th highest rank of a positive item, and
    // cPos is the total number of positive items.

    const int iRankItemPos  = ranker.GetRank(iItemPos);
    const int iRankItemNeg  = ranker.GetRank(iItemNeg);

    // Search for the position of the two items to swap
    const vector<int>::iterator itItemPos = upper_bound(rankpos_vec_.begin(), rankpos_vec_.begin() + cPos, iRankItemPos);
    const vector<int>::iterator itItemNeg = upper_bound(rankpos_vec_.begin(), rankpos_vec_.begin() + cPos, iRankItemNeg);

    // The number of positive items up to and including iItemPos
    const int cNumPosNotBelowItemPos = (int)(itItemPos - rankpos_vec_.begin());

    // The number of positive items up to iItemNeg (Note: Cannot include iItemNeg itself)
    const unsigned int cNumPosAboveItemNeg    = (unsigned int)(itItemNeg  - rankpos_vec_.begin());

    // Range of indices of positive items between iRankItemPos and iRankItemNeg (exclusively)
    int cIntermediateHigh, cIntermediateLow;

    // Current contribution of iItemPos
    double dContribBefore =  (double) cNumPosNotBelowItemPos / iRankItemPos;

    double dSign, dContribAfter;

    if (iRankItemNeg > iRankItemPos)
    {
        // MAP is decreasing
        dSign = -1.0;

        // The first positive item after iRankItemPos
        cIntermediateLow = cNumPosNotBelowItemPos;

        // The last positive item before iRankItemNeg
        cIntermediateHigh = cNumPosAboveItemNeg - 1;

        // Note: iItemPos already counted in cNumPosAboveItemNeg
        dContribAfter = (double)cNumPosAboveItemNeg / iRankItemNeg;
    }
    else
    {
        // MAP is increasing
        dSign = 1.0;

        // The first positive result after iRankItemNeg
        cIntermediateLow = cNumPosAboveItemNeg;

        // The first positive result after iRankItemPos, minus iItemPos itself
        cIntermediateHigh = cNumPosNotBelowItemPos - 2;

        // Note: iItemPos not yet counted in cNumPosAboveItemNeg
        dContribAfter = (double) (cNumPosAboveItemNeg + 1) / iRankItemNeg;
    }

    // The direct effect of switching iItemPos
    double dDiff = dContribAfter - dContribBefore;

    // The indirect effect for all items in between the two items
    for (int j = cIntermediateLow; j <= cIntermediateHigh; j++)
    {
        dDiff += dSign / rankpos_vec_[j];
    }

    return dDiff / cPos;
}


double CMAP::Measure(const double* const adY, const CRanker& ranker)
{
    unsigned int cPos;

    SortRankPos(adY, ranker, rankpos_vec_, cPos);

    if (cPos == 0)
    {
        return 0.0;
    }

    // Now veccRankPos[i] is the i-th highest rank of a positive item

    double dPrec = 0.0;
    for (unsigned int j = 0; j < cPos; j++)
    {
        dPrec += double(j + 1) / rankpos_vec_[j];
    }

    return dPrec / cPos;
}


CPairwise::CPairwise(const double* adgroups, const char* szIRMeasure, int cTrain)
{

	// Set up adGroup - this is not required
	kGroups_ = adgroups;

	// Set up the number of groups - this used externally
	SetNumGroups(GBM_FUNC::NumGroups(adgroups, cTrain));

    // Construct the IR Measure
    if (!strcmp(szIRMeasure, "conc"))
    {
      pirm_.reset(new CConc());
    }
    else if (!strcmp(szIRMeasure, "map"))
    {
      pirm_.reset(new CMAP());
    }
    else if (!strcmp(szIRMeasure, "mrr"))
    {
      pirm_.reset(new CMRR());
    }
    else
      {
        if (strcmp(szIRMeasure, "ndcg"))
        {
            Rprintf("Unknown IR measure '%s' in initialization, using 'ndcg' instead\n", szIRMeasure);
        }
        pirm_.reset(new CNDCG());
      }
}

CDistribution* CPairwise::Create(DataDistParams& distParams)
{
	// Create pointers to pairwise
	Rcpp::NumericVector miscVec(distParams.misc[0]);
	const double* adgroup = 0;

	if(!GBM_FUNC::has_value(miscVec))
	{
		throw GBM::Failure("Pairwise requires misc to initialize");
	}
	else
	{
		adgroup = miscVec.begin();
	}
	return new CPairwise(adgroup, distParams.irmeasure, distParams.num_trainrows);
}

CPairwise::~CPairwise()
{
}


// Auxiliary function for addition of optional offset parameter
inline const double* OffsetVector(const double* const adX, const double* const adOffset, unsigned int iStart, unsigned int iEnd, vector<double>& vecBuffer)
{
    if (adOffset == NULL)
    {
        // Optional second argument is not set, just return first one
        return adX + iStart;
    }
    else
    {
        for (unsigned int i = iStart, iOut = 0; i < iEnd; i++, iOut++)
        {
            vecBuffer[iOut] = adX[i] + adOffset[i];
        }
        return &vecBuffer[0];
    }
}

void CPairwise::ComputeWorkingResponse
(
 const CDataset& data,
 const double *adF,
 double *adZ
)
{
#ifdef NOISY_DEBUG
    Rprintf("compute working response, nTrain = %u\n", nTrain);
#endif
    
    if (data.get_trainsize() <= 0) return;
    
    // Iterate through all groups, compute gradients
    
    unsigned int iItemStart = 0;
    unsigned int iItemEnd   = 0;
    
    while (iItemStart < data.get_trainsize())
      {
	adZ[iItemEnd]           = 0;
	hessian_[iItemEnd]   = 0;
	
	const double dGroup = kGroups_[iItemStart];
	
	// Find end of current group, initialize working response
	for (iItemEnd = iItemStart + 1; iItemEnd < data.get_trainsize() && kGroups_[iItemEnd] == dGroup; iItemEnd++)
	  {
	    // Clear gradients from last iteration
	    adZ[iItemEnd]         = 0;
	    hessian_[iItemEnd] = 0;
	  }
	
#ifdef NOISY_DEBUG
	// Check sorting
	for (unsigned int i = iItemStart; i < iItemEnd-1; i++) {
	  if (data. y_ptr()[i] < data. y_ptr()[i+1]) {
	    throw GBM::Failure("sorting failed in pairwise?");
	  }
	}  
#endif

	if (data.get_bag_element(iItemStart))
	  {
	    // Group is part of the training set
	    
	    const int cNumItems = iItemEnd - iItemStart;
	    
	    // If offset given, add up current scores
	    const double* adFPlusOffset = OffsetVector(adF, data.offset_ptr(), iItemStart, iItemEnd, func_est_plus_offset_);
	    
	    // Accumulate gradients
	    // TODO: Implement better way to ensure casting robust to overflow
	    int intGroup = 0;
	    if(fabs(dGroup) > nextafter(INT_MAX, 0) || isnan(dGroup))
	    {
	    	intGroup = copysign(INT_MAX, dGroup);
	    }
	    else
	    {
	    	intGroup = (int) dGroup;
	    }

	    ComputeLambdas(intGroup, cNumItems, data.y_ptr() + iItemStart, adFPlusOffset, data.weight_ptr() + iItemStart, adZ + iItemStart, &hessian_[iItemStart]);
	  }
	
	// Next group
	iItemStart = iItemEnd;
      }
}


// Referring to MSR-TR-2010-82-2, section 7 (see also the vignette):
//
// Let P be the set of pairs (i,j) where Y(i)>Y(j) (i is better than j).
// The approximation to the IR measure is the utility function C (to be maximized)
//   C
//   = \Sum_{(i,j) in P} |Delta Z_ij| C(s_i - s_j)
//   = \Sum_{(i,j) in P} |Delta Z_ij| / (1 + std::exp(-(s_i - s_j))),
// where |Delta Z_ij| is the cost of swapping (only) i and j in the current ranking,
// and s_i, s_j are the prediction scores (sum of the tree predictions) for items
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

void CPairwise::ComputeLambdas(int iGroup, unsigned int cNumItems, const double* const adY, const double* const adF, const double* const adWeight, double* adZ, double* adDeriv)
{
    // Assumption: Weights are constant within group
    if (adWeight[0] <= 0)
    {
        return;
    }

    // Normalize for maximum achievable group score
    const double dMaxScore = pirm_->MaxMeasure(iGroup, adY, cNumItems);

    if (dMaxScore <= 0.0)
    {
        // No pairs
        return;
    }

    // Rank items by current score
    ranker_.SetGroupScores(adF, cNumItems);
    ranker_.Rank();

    double dLabelCurrent    = adY[0];

    // First index of instance that has dLabelCurrent
    // (i.e., each smaller index corresponds to better item)
    unsigned int iLabelCurrentStart  = 0;

    // Number of pairs with unequal labels
    unsigned int cPairs   = 0;

#ifdef NOISY_DEBUG
    double dMeasureBefore = pirm_->Measure(adY, ranker_);
#endif

    for (unsigned int j = 1; j < cNumItems; j++)
    {
        const double dYj     = adY[j];

        if (dYj != dLabelCurrent)
        {
            iLabelCurrentStart = j;
            dLabelCurrent      = dYj;
        }

        for (unsigned int i = 0; i < iLabelCurrentStart; i++)
        {
            // Instance i is better than j

            const double dSwapCost = fabs(pirm_->SwapCost(i, j, adY, ranker_));

#ifdef NOISY_DEBUG
            double dDelta    = fabs(pirm_->SwapCost(i, j, adY, ranker_));
            const int cRanki = ranker_.GetRank(i);
            const int cRankj = ranker_.GetRank(j);
            ranker_.SetRank(i, cRankj);
            ranker_.SetRank(j, cRanki);
            double dMeasureAfter = pirm_->Measure(adY, ranker_);

            if (fabs(dMeasureBefore-dMeasureAfter) - dDelta > 1e-5)
            {
                Rprintf("%f %f %f %f %f %d %d\n", pirm_->SwapCost(i, j, adY, ranker_), dMeasureBefore, dMeasureAfter, dMeasureBefore - dMeasureAfter, dDelta , i, j);
                for (unsigned int k = 0; k < num_items_; k++)
                {
                    Rprintf("%d\t%d\t%f\t%f\n", k, ranker_.GetRank(k), adY[k], adF[k]);
                }
		throw GBM::Failure("the impossible happened");
            }
	    if (fabs(dMeasureBefore - dMeasureAfter) - fabs(dDelta) >= 1e-5) {
	      throw GBM::Failure("the impossible happened");
	    }
            ranker_.SetRank(j, cRankj);
            ranker_.SetRank(i, cRanki);
	    
#endif
	    if (!isfinite(dSwapCost)) {
	      throw GBM::Failure("infinite swap cost");
	    }

            if (dSwapCost > 0.0)
            {
                cPairs++;
                const double dRhoij    = 1.0 / (1.0 + std::exp(adF[i]- adF[j])) ;
                if (!isfinite(dRhoij)) {
		  throw GBM::Failure("unanticipated infinity");
		};

                const double dLambdaij = dSwapCost * dRhoij;
                adZ[i] += dLambdaij;
                adZ[j] -= dLambdaij;
                const double dDerivij  = dLambdaij * (1.0 - dRhoij);
		if (dDerivij < 0) {
		  throw GBM::Failure("negative derivative!");
		}
                adDeriv[i] += dDerivij;
                adDeriv[j] += dDerivij;
            }
        }
    }

    if (cPairs > 0)
    {
        // Normalize for number of training pairs
        const double dQNorm     = 1.0 / (dMaxScore * cPairs);

        for (unsigned int j = 0; j < cNumItems; j++)
        {
            adZ[j]     *= dQNorm;
            adDeriv[j] *= dQNorm;
        }
    }
}

void CPairwise::Initialize
(
	const CDataset& data
)
{
  if (data.nrow() <= 0) return;
  
  // Allocate memory for derivative buffer
  hessian_.resize(data.nrow());

  // Count the groups and number of items per group
  unsigned int cMaxItemsPerGroup = 0;
  double       dMaxGroup         = 0;
  
  unsigned int iItemStart        = 0;
  unsigned int iItemEnd          = 0;
  
  while (iItemStart < data.nrow())
    {
      
      const double dGroup = kGroups_[iItemStart];
      
      // Find end of current group
      for (iItemEnd = iItemStart + 1; iItemEnd < data.nrow() && kGroups_[iItemEnd] == dGroup; iItemEnd++);
      
      const unsigned int cNumItems = iItemEnd - iItemStart;
      if (cNumItems > cMaxItemsPerGroup)
	{
	  cMaxItemsPerGroup = cNumItems;
	}
      if (dGroup > dMaxGroup)
	{
	  dMaxGroup = dGroup;
	}
      
      // Next group
      iItemStart = iItemEnd;
    }
  
  // Allocate buffer for offset addition
  func_est_plus_offset_.resize(cMaxItemsPerGroup);
  
  // Allocate ranker memory
  ranker_.Init(cMaxItemsPerGroup);
  
  // Allocate IR measure memory
  
  // The last element of adGroup specifies the cutoff
  // (zero means no cutoff)
  unsigned int cRankCutoff = cMaxItemsPerGroup;
  if (kGroups_[data.nrow()] > 0)
    {
      cRankCutoff = (unsigned int)kGroups_[data.nrow()];
    }

  	// TODO: Make More robust against overflow
  	unsigned long ulMaxGroup = 0;
	if(fabs(dMaxGroup) > nextafter(ULONG_MAX, 0) || isnan(dMaxGroup))
	{
		ulMaxGroup = copysign(ULONG_MAX, dMaxGroup);
	}
	else
	{
		ulMaxGroup = (unsigned long) dMaxGroup;
	}
  pirm_->Init(ulMaxGroup, cMaxItemsPerGroup, cRankCutoff);
#ifdef NOISY_DEBUG
  Rprintf("Initialization: instances=%ld, groups=%u, max items per group=%u, rank cutoff=%u, offset specified: %d\n", cLength, (unsigned long)dMaxGroup, cMaxItemsPerGroup, cRankCutoff, (data.offset_ptr() != NULL));
#endif
}

double CPairwise::InitF
(
	const CDataset& data
)
{
    return 0.0;
}


double CPairwise::Deviance
(
   const CDataset& data,
   const double *adF,
   bool isValidationSet
)
{

    // Shift adGroup to validation set if necessary
	long cLength = data.get_trainsize();
    if(isValidationSet)
    {
    	cLength = data.get_validsize();
    	data.shift_to_validation();
    	kGroups_=shift_ptr(kGroups_, data.get_trainsize());

    }

    if (cLength <= 0)
	{
    	// NB: SWITCH BACK TO TRAIN BEFORE LEAVING
    	data.shift_to_train();
    	kGroups_=shift_ptr(kGroups_, -(data.get_trainsize()));
    	return 0;
	}

    double dL = 0.0;
    double dW = 0.0;

    unsigned int iItemStart  = 0;
    unsigned int iItemEnd    = iItemStart;
    const unsigned int cEnd = cLength;

    while (iItemStart < cEnd)
    {
        const double dGroup = kGroups_[iItemStart];
        const double dWi    = data.weight_ptr()[iItemStart];

        // Find end of current group
        for (iItemEnd = iItemStart + 1; iItemEnd < cEnd && kGroups_[iItemEnd] == dGroup; iItemEnd++) ;

        const int cNumItems = iItemEnd - iItemStart;
        // TODO: Implement better way to ensure casting robust to overflow
		int intGroup = 0;
		if(fabs(dGroup) > nextafter(INT_MAX, 0) || isnan(dGroup))
		{
			intGroup = copysign(INT_MAX, dGroup);
		}
		else
		{
			intGroup = (int) dGroup;
		}
        const double dMaxScore = pirm_->MaxMeasure(intGroup, data.y_ptr() + iItemStart, cNumItems);

        if (dMaxScore > 0.0)
        {
            // Rank items by current score

            // If offset given, add up current scores
            const double* adFPlusOffset = OffsetVector(adF, data.offset_ptr(), iItemStart, iItemEnd, func_est_plus_offset_);

            ranker_.SetGroupScores(adFPlusOffset, cNumItems);
            ranker_.Rank();

            dL += dWi * pirm_->Measure(data.y_ptr() + iItemStart, ranker_) / dMaxScore;
            dW += dWi;
        }
        // Next group
        iItemStart = iItemEnd;
    }

    // Reset adGroup if required
    if(isValidationSet)
    {
    	data.shift_to_train();
    	kGroups_=shift_ptr(kGroups_, -(data.get_trainsize()));

    }

   // Loss = 1 - utility
   return 1.0 - dL / dW;
}


void CPairwise::FitBestConstant
(
	const CDataset& data,
    const double *adF,
    unsigned long cTermNodes,
    double* adZ,
    CTreeComps& treeComps
)
{

#ifdef NOISY_DEBUG
    Rprintf("FitBestConstant, nTrain = %u,  cTermNodes = %d, \n", nTrain, cTermNodes);
#endif

    // Assumption: ComputeWorkingResponse() has been executed before with
    // the same arguments

    // Allocate space for numerators and denominators, and set to zero
    fit_numerator_.reserve(cTermNodes);
    fit_denominator_.reserve(cTermNodes);
    for (unsigned int i = 0; i < cTermNodes; i++)
      {
	fit_numerator_[i]   = 0.0;
	fit_denominator_[i] = 0.0;
      }

    for (unsigned int iObs = 0; iObs < data.get_trainsize(); iObs++)
    {
      if (data.get_bag_element(iObs))
        {
#ifdef NOISY_DEBUG
          if (!(isfinite(data.weight_ptr()[iObs]) &&
		isfinite(adZ[iObs]) &&
		isfinite(hessian_[iObs]))) {
	    throw GBM::Failure("unanticipated infinities");
	  };
#endif

            fit_numerator_[treeComps.get_node_assignments()[iObs]]   += data.weight_ptr()[iObs] * adZ[iObs];
            fit_denominator_[treeComps.get_node_assignments()[iObs]] += data.weight_ptr()[iObs] * hessian_[iObs];
        }
    }

    for (unsigned int iNode = 0; iNode < cTermNodes; iNode++)
    {
        if (treeComps.get_terminal_nodes()[iNode] != NULL)
        {
            if (fit_denominator_[iNode] <= 0.0)
            {
            	treeComps.get_terminal_nodes()[iNode]->prediction = 0.0;
            }
            else
            {
            	treeComps.get_terminal_nodes()[iNode]->prediction =
                    fit_numerator_[iNode]/fit_denominator_[iNode];
            }
        }
    }
}


double CPairwise::BagImprovement
(
	const CDataset& data,
    const double *adF,
    const double shrinkage,
    const double* adFadj
)
{

#ifdef NOISY_DEBUG
    Rprintf("BagImprovement, nTrain = %u\n", nTrain);
#endif

    if (data.get_trainsize() <= 0)
    {
        return 0;
    }

    double dL = 0.0;
    double dW = 0.0;

    unsigned int iItemStart = 0;
    unsigned int iItemEnd   = 0;


    while (iItemStart < data.get_trainsize())
    {
        const double dGroup = kGroups_[iItemStart];

        // Find end of current group
        for (iItemEnd = iItemStart + 1; iItemEnd < data.get_trainsize() && kGroups_[iItemEnd] == dGroup; iItemEnd++) ;

        if (!data.get_bag_element(iItemStart))
        {
            // Group was held out of training set

            const unsigned int cNumItems = iItemEnd - iItemStart;
            // TODO: Implement better way to ensure casting robust to overflow
			int intGroup = 0;
			if(fabs(dGroup) > nextafter(INT_MAX, 0) || isnan(dGroup))
			{
				intGroup = copysign(INT_MAX, dGroup);
			}
			else
			{
				intGroup = (int) dGroup;
			}
            const double dMaxScore = pirm_->MaxMeasure(intGroup, data.y_ptr() + iItemStart, cNumItems);

            if (dMaxScore > 0.0)
            {
                // If offset given, add up current scores
                const double* adFPlusOffset = OffsetVector(adF, data.offset_ptr(), iItemStart, iItemEnd, func_est_plus_offset_);

                // Compute score according to old score, adF
                ranker_.SetGroupScores(adFPlusOffset, cNumItems);
                ranker_.Rank();
                const double dOldScore = pirm_->Measure(data.y_ptr() + iItemStart, ranker_);

                // Compute score according to new score: adF' =  adF + dStepSize * adFadj
                for (unsigned int i = 0; i < cNumItems; i++)
                {
                    ranker_.AddToScore(i, adFadj[i+iItemStart] * shrinkage);
                }

                const double dWi = data.weight_ptr()[iItemStart];

                if (ranker_.Rank())
                {
                    // Ranking changed
                    const double dNewScore = pirm_->Measure(data.y_ptr() + iItemStart, ranker_);
                    dL                    += dWi * (dNewScore - dOldScore) / dMaxScore;
                }
                dW += dWi;
            }
        }

        // Next group
        iItemStart = iItemEnd;

    }

    return dL / dW;

}

void CPairwise::BagData(CDataset& data) {
  double dLastGroup = -1;
  bool fChosen = false;
  unsigned int cBagged = 0;
  unsigned int cBaggedGroups = 0;
  unsigned int cSeenGroups   = 0;
  unsigned int cTotalGroupsInBag = (unsigned long)(data.get_bagfraction() * GetNumGroups());

  if (cTotalGroupsInBag <= 0) {
    cTotalGroupsInBag = 1;
  }
  
  for(unsigned long i=0; i< data.get_trainsize(); i++) {
    const double dGroup = kGroups_[i];
    
    if(dGroup != dLastGroup) {
      if (cBaggedGroups >= cTotalGroupsInBag) {
	break;
      }
      
      // Group changed, make a new decision
      fChosen = (unif_rand()*(GetNumGroups() - cSeenGroups) <
		 cTotalGroupsInBag - cBaggedGroups);
      if(fChosen) {
	cBaggedGroups++;
      }
      dLastGroup = dGroup;
      cSeenGroups++;
    }
    
    if(fChosen) {
      data.set_bag_element(i);
      cBagged++;
    }
  }
  
}
