// Implementation file for 'pairwise' distribution
//
// Author: Stefan Schroedl (schroedl@a9.com)

#include "pairwise.h"
#include "gbmFunc.h"
#include <iostream>
#include <vector>
#include <algorithm>

//#define NOISY_DEBUG
#ifdef NOISY_DEBUG
#endif

void CRanker::Init(unsigned int cMaxItemsPerGroup)
{
    // Allocate sorting buffers
    vecdipScoreRank.resize(cMaxItemsPerGroup);
    vecpdipScoreRank.resize(cMaxItemsPerGroup);
}

bool CRanker::SetGroupScores(const double* const adScores, const unsigned int cNumItems)
{
    const double dEPS = 1e-10;

    if (cNumItems > vecdipScoreRank.size())
    {
        // Allocate additional space
        // (We should never get here if CPairwise::Initialize has been called before, as expected)
        Init(cNumItems);
    }
    this->cNumItems = cNumItems;

    // Copy scores to buffer, and
    // initialize pointer array to score entries

    for(unsigned int i = 0; i < cNumItems; i++)
    {
        // Add small random number to break possible ties
        vecdipScoreRank[i].first = adScores[i] + dEPS * (unif_rand() - 0.5);

        vecpdipScoreRank[i] = &(vecdipScoreRank[i]);
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

    sort(vecpdipScoreRank.begin(), vecpdipScoreRank.begin() + cNumItems, comp);

    bool bChanged = false;

    // Create inverted rank lookup

    for(unsigned int i = 0; i < cNumItems; i++)
    {
        // Note: ranks are 1-based
        const unsigned int cNewRank = i + 1;
        if (!bChanged)
        {
            bChanged = (cNewRank != vecpdipScoreRank[i]->second);
        }
        // Store the rank with the corresponding score in the vecdipScoreRank array
        vecpdipScoreRank[i]->second = cNewRank;
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
    veccPairCount.resize(cMaxGroup + 1, -1);
}

unsigned int CConc::PairCount(unsigned int iGroup, const double* const adY, unsigned int cNumItems)
{
    if (iGroup >= veccPairCount.size())
    {
        // Allocate additional space
        // (We should never get here if CPairwise::Initialize has been called before, as expected)
        veccPairCount.resize(iGroup + 1, -1);
    }

    if (veccPairCount[iGroup] < 0.0)
    {
        // Not yet initialized
        veccPairCount[iGroup] = ComputePairCount(adY, cNumItems);
    }
    return veccPairCount[iGroup];
}

// Calculate the number of pairs with different labels, and store in veccPairCount
// Assumption: instances are sorted such that labels are non-increasing
int CConc::ComputePairCount(const double* const adY, unsigned int cNumItems)
{
    if (!AnyPairs(adY, cNumItems))
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

    vecdRankWeight.resize(cMaxItemsPerGroup + 1, 0.0);

    const unsigned int cMaxRank = std::min((unsigned int)cMaxItemsPerGroup, GetCutoffRank());

    // Precompute rank weights
    for (unsigned int i = 1; i <= cMaxRank; i++)
    {
        vecdRankWeight[i] = std::log((double)2) / log((double)(i+1));
    }

    // Allocate buffer
    vecdMaxDCG.resize(cMaxGroup + 1, -1.0);
}


// Sum of target values, weighted by rank weight
double CNDCG::Measure(const double* const adY, const CRanker& ranker)
{
    double dScore = 0;

    for (unsigned int i = 0; i < ranker.GetNumItems(); i++)
    {
        dScore += adY[i] * vecdRankWeight[ranker.GetRank(i)];
    }
    return dScore;
}

double CNDCG::MaxMeasure(unsigned int iGroup, const double* const adY, unsigned int cNumItems)
{
    if (iGroup >= vecdMaxDCG.size())
    {
        // Allocate additional space
        // (We should never get here if CPairwise::Initialize has been called before, as expected)
    	// JHickey: Then should we call initialize in the constructor?
        vecdMaxDCG.resize(iGroup + 1, -1.0);
    }

    if (vecdMaxDCG[iGroup] < 0.0)
    {
        // Not initialized

        if (!AnyPairs(adY, cNumItems))
        {
            // No training pairs exist
            vecdMaxDCG[iGroup] = 0.0;
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
                dScore += adY[i] * vecdRankWeight[i + 1];
                i++;
            }

            vecdMaxDCG[iGroup] = dScore;

#ifdef NOISY_DEBUG
            if (vecdMaxDCG[iGroup] == 0)
            {
                Rprintf("max score is 0: iGroup = %d, maxScore = %f\n", 
                        iGroup,  vecdMaxDCG[iGroup]);
                throw GBM::failure();
            }
#endif
        }
    }

    return vecdMaxDCG[iGroup];
}


double CNDCG::SwapCost(int iItemBetter, int iItemWorse, const double* const adY, const CRanker& ranker) const
{
    const unsigned int cRanki = ranker.GetRank(iItemBetter);
    const unsigned int cRankj = ranker.GetRank(iItemWorse);
    return (vecdRankWeight[cRanki] - vecdRankWeight[cRankj]) * (adY[iItemBetter] - adY[iItemWorse]);
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

    const unsigned int cNumItems =  min(ranker.GetNumItems(), GetCutoffRank());

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

    const unsigned int cCutoffRank = GetCutoffRank();
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
    veccRankPos.resize(cMaxItemsPerGroup + 1);
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

    SortRankPos(adY, ranker, veccRankPos, cPos);

    if (cPos == 0)
    {
        return 0.0;
    }

    // Now veccRankPos[i] is the i-th highest rank of a positive item, and
    // cPos is the total number of positive items.

    const int iRankItemPos  = ranker.GetRank(iItemPos);
    const int iRankItemNeg  = ranker.GetRank(iItemNeg);

    // Search for the position of the two items to swap
    const vector<int>::iterator itItemPos = upper_bound(veccRankPos.begin(), veccRankPos.begin() + cPos, iRankItemPos);
    const vector<int>::iterator itItemNeg = upper_bound(veccRankPos.begin(), veccRankPos.begin() + cPos, iRankItemNeg);

    // The number of positive items up to and including iItemPos
    const unsigned int cNumPosNotBelowItemPos = (unsigned int)(itItemPos - veccRankPos.begin());

    // The number of positive items up to iItemNeg (Note: Cannot include iItemNeg itself)
    const unsigned int cNumPosAboveItemNeg    = (unsigned int)(itItemNeg  - veccRankPos.begin());

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
        dDiff += dSign / veccRankPos[j];
    }

    return dDiff / cPos;
}


double CMAP::Measure(const double* const adY, const CRanker& ranker)
{
    unsigned int cPos;

    SortRankPos(adY, ranker, veccRankPos, cPos);

    if (cPos == 0)
    {
        return 0.0;
    }

    // Now veccRankPos[i] is the i-th highest rank of a positive item

    double dPrec = 0.0;
    for (unsigned int j = 0; j < cPos; j++)
    {
        dPrec += double(j + 1) / veccRankPos[j];
    }

    return dPrec / cPos;
}


CPairwise::CPairwise(SEXP radMisc,
					const char* szIRMeasure, int& cTrain): CDistribution(radMisc)
{

	// Set up adGroup - this is not required
	adGroup = CDistribution::misc_ptr(false);

	// Set up the number of groups - this used externally
	SetNumGroups(GBM_FUNC::numGroups(CDistribution::misc_ptr(true), cTrain));

    // Construct the IR Measure
    if (!strcmp(szIRMeasure, "conc"))
    {
      pirm.reset(new CConc());
    }
    else if (!strcmp(szIRMeasure, "map"))
    {
      pirm.reset(new CMAP());
    }
    else if (!strcmp(szIRMeasure, "mrr"))
    {
      pirm.reset(new CMRR());
    }
    else
      {
        if (strcmp(szIRMeasure, "ndcg"))
        {
            Rprintf("Unknown IR measure '%s' in initialization, using 'ndcg' instead\n", szIRMeasure);
        }
        pirm.reset(new CNDCG());
      }
}

CDistribution* CPairwise::Create(SEXP radMisc,
										const char* szIRMeasure, int& cTrain)
{

	return new CPairwise(radMisc, szIRMeasure, cTrain);
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
 const CDataset* pData,
 const double *adF,
 double *adZ
)
{
#ifdef NOISY_DEBUG
    Rprintf("compute working response, nTrain = %u\n", nTrain);
#endif
    
    if (pData->get_trainSize() <= 0) return;
    
    // Iterate through all groups, compute gradients
    
    unsigned int iItemStart = 0;
    unsigned int iItemEnd   = 0;
    
    while (iItemStart < pData->get_trainSize())
      {
	adZ[iItemEnd]           = 0;
	vecdHessian[iItemEnd]   = 0;
	
	const double dGroup = adGroup[iItemStart];
	
	// Find end of current group, initialize working response
	for (iItemEnd = iItemStart + 1; iItemEnd < pData->get_trainSize() && adGroup[iItemEnd] == dGroup; iItemEnd++)
	  {
	    // Clear gradients from last iteration
	    adZ[iItemEnd]         = 0;
	    vecdHessian[iItemEnd] = 0;
	  }
	
#ifdef NOISY_DEBUG
	// Check sorting
	for (unsigned int i = iItemStart; i < iItemEnd-1; i++) {
	  if (pData-> y_ptr()[i] < pData-> y_ptr()[i+1]) {
	    throw GBM::failure("sorting failed in pairwise?");
	  }
	}  
#endif

	if (pData->GetBag()[iItemStart])
	  {
	    // Group is part of the training set
	    
	    const int cNumItems = iItemEnd - iItemStart;
	    
	    // If offset given, add up current scores
	    const double* adFPlusOffset = OffsetVector(adF, pData->offset_ptr(false), iItemStart, iItemEnd, vecdFPlusOffset);
	    
	    // Accumulate gradients
	    ComputeLambdas((int)dGroup, cNumItems, pData->y_ptr() + iItemStart, adFPlusOffset, pData->weight_ptr() + iItemStart, adZ + iItemStart, &vecdHessian[iItemStart]);
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
    const double dMaxScore = pirm->MaxMeasure(iGroup, adY, cNumItems);

    if (dMaxScore <= 0.0)
    {
        // No pairs
        return;
    }

    // Rank items by current score
    ranker.SetGroupScores(adF, cNumItems);
    ranker.Rank();

    double dLabelCurrent    = adY[0];

    // First index of instance that has dLabelCurrent
    // (i.e., each smaller index corresponds to better item)
    unsigned int iLabelCurrentStart  = 0;

    // Number of pairs with unequal labels
    unsigned int cPairs   = 0;

#ifdef NOISY_DEBUG
    double dMeasureBefore = pirm->Measure(adY, ranker);
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

            const double dSwapCost = fabs(pirm->SwapCost(i, j, adY, ranker));

#ifdef NOISY_DEBUG
            double dDelta    = fabs(pirm->SwapCost(i, j, adY, ranker));
            const int cRanki = ranker.GetRank(i);
            const int cRankj = ranker.GetRank(j);
            ranker.SetRank(i, cRankj);
            ranker.SetRank(j, cRanki);
            double dMeasureAfter = pirm->Measure(adY, ranker);

            if (fabs(dMeasureBefore-dMeasureAfter) - dDelta > 1e-5)
            {
                Rprintf("%f %f %f %f %f %d %d\n", pirm->SwapCost(i, j, adY, ranker), dMeasureBefore, dMeasureAfter, dMeasureBefore - dMeasureAfter, dDelta , i, j);
                for (unsigned int k = 0; k < cNumItems; k++)
                {
                    Rprintf("%d\t%d\t%f\t%f\n", k, ranker.GetRank(k), adY[k], adF[k]);
                }
		throw GBM::failure("the impossible happened");
            }
	    if (fabs(dMeasureBefore - dMeasureAfter) - fabs(dDelta) >= 1e-5) {
	      throw GBM::failure("the impossible happened");
	    }
            ranker.SetRank(j, cRankj);
            ranker.SetRank(i, cRanki);
	    
	    if (!isfinite(dSwapCost)) {
	      throw GBM::failure("infinite swap cost");
	    }
#endif

            if (dSwapCost > 0.0)
            {
#ifdef NOISY_DEBUG
                cPairs++;
                const double dRhoij    = 1.0 / (1.0 + std::exp(adF[i]- adF[j])) ;
                if (!isfinite(dRhoij)) {
		  throw GBM::failure("unanticipated infinity");
		};

                const double dLambdaij = dSwapCost * dRhoij;
                adZ[i] += dLambdaij;
                adZ[j] -= dLambdaij;
                const double dDerivij  = dLambdaij * (1.0 - dRhoij);
		if (dDerivij < 0) {
		  throw GBM::failure("negative derivative!");
		}
                adDeriv[i] += dDerivij;
                adDeriv[j] += dDerivij;
#endif
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
	const CDataset* pData
)
{
  if (pData->nrow() <= 0) return;
  
  // Allocate memory for derivative buffer
  vecdHessian.resize(pData->nrow());

  // Count the groups and number of items per group
  unsigned int cMaxItemsPerGroup = 0;
  double       dMaxGroup         = 0;
  
  unsigned int iItemStart        = 0;
  unsigned int iItemEnd          = 0;
  
  while (iItemStart < pData->nrow())
    {
      
      const double dGroup = adGroup[iItemStart];
      
      // Find end of current group
      for (iItemEnd = iItemStart + 1; iItemEnd < pData->nrow() && adGroup[iItemEnd] == dGroup; iItemEnd++);
      
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
  vecdFPlusOffset.resize(cMaxItemsPerGroup);
  
  // Allocate ranker memory
  ranker.Init(cMaxItemsPerGroup);
  
  // Allocate IR measure memory
  
  // The last element of adGroup specifies the cutoff
  // (zero means no cutoff)
  unsigned int cRankCutoff = cMaxItemsPerGroup;
  if (adGroup[pData->nrow()] > 0)
    {
      cRankCutoff = (unsigned int)adGroup[pData->nrow()];
    }
  pirm->Init((unsigned long)dMaxGroup, cMaxItemsPerGroup, cRankCutoff);
#ifdef NOISY_DEBUG
  Rprintf("Initialization: instances=%ld, groups=%u, max items per group=%u, rank cutoff=%u, offset specified: %d\n", cLength, (unsigned long)dMaxGroup, cMaxItemsPerGroup, cRankCutoff, (pData->offset_ptr(false) != NULL));
#endif
}

double CPairwise::InitF
(
	const CDataset* pData
)
{
    return 0.0;
}


double CPairwise::Deviance
(
   const CDataset* pData,
   const double *adF,
   bool isValidationSet
)
{

    // Shift adGroup to validation set if necessary
	long cLength = pData->get_trainSize();
	if (cLength <= 0)
	{
	        return 0;
	}

    if(isValidationSet)
    {
    	pData->shift_to_validation();
    	adGroup=shift_ptr(CDistribution::misc_ptr(false), pData->get_trainSize());
    	cLength = pData->GetValidSize();
    }


    double dL = 0.0;
    double dW = 0.0;

    unsigned int iItemStart  = 0;
    unsigned int iItemEnd    = iItemStart;
    const unsigned int cEnd = cLength;




    while (iItemStart < cEnd)
    {
        const double dGroup = adGroup[iItemStart];
        const double dWi    = pData->weight_ptr()[iItemStart];

        // Find end of current group
        for (iItemEnd = iItemStart + 1; iItemEnd < cEnd && adGroup[iItemEnd] == dGroup; iItemEnd++) ;

        const int cNumItems = iItemEnd - iItemStart;

        const double dMaxScore = pirm->MaxMeasure((int)dGroup, pData->y_ptr() + iItemStart, cNumItems);

        if (dMaxScore > 0.0)
        {
            // Rank items by current score

            // If offset given, add up current scores
            const double* adFPlusOffset = OffsetVector(adF, pData->offset_ptr(false), iItemStart, iItemEnd, vecdFPlusOffset);

            ranker.SetGroupScores(adFPlusOffset, cNumItems);
            ranker.Rank();

            dL += dWi * pirm->Measure(pData->y_ptr() + iItemStart, ranker) / dMaxScore;
            dW += dWi;
        }
        // Next group
        iItemStart = iItemEnd;
    }

    // Reset adGroup if required
    if(isValidationSet)
    {
    	pData->shift_to_train();
    	adGroup = CDistribution::misc_ptr(false);
    }
   // Loss = 1 - utility
   return 1.0 - dL / dW;
}


void CPairwise::FitBestConstant
(
	const CDataset* pData,
    const double *adF,
    unsigned long cTermNodes,
    double* adZ,
    CTreeComps* pTreeComps
)
{

#ifdef NOISY_DEBUG
    Rprintf("FitBestConstant, nTrain = %u,  cTermNodes = %d, \n", nTrain, cTermNodes);
#endif

    // Assumption: ComputeWorkingResponse() has been executed before with
    // the same arguments

    // Allocate space for numerators and denominators, and set to zero
    vecdNum.reserve(cTermNodes);
    vecdDenom.reserve(cTermNodes);
    for (unsigned int i = 0; i < cTermNodes; i++)
      {
	vecdNum[i]   = 0.0;
	vecdDenom[i] = 0.0;
      }

    for (unsigned int iObs = 0; iObs < pData->get_trainSize(); iObs++)
    {
      if (pData->GetBagElem(iObs))
        {
#ifdef NOISY_DEBUG
          if (!(isfinite(pData->weight_ptr()[iObs]) &&
		isfinite(adZ[iObs]) &&
		isfinite(vecdHessian[iObs]))) {
	    throw GBM::failure("unanticipated infinities");
	  };
#endif

            vecdNum[pTreeComps->GetNodeAssign()[iObs]]   += pData->weight_ptr()[iObs] * adZ[iObs];
            vecdDenom[pTreeComps->GetNodeAssign()[iObs]] += pData->weight_ptr()[iObs] * vecdHessian[iObs];
        }
    }

    for (unsigned int iNode = 0; iNode < cTermNodes; iNode++)
    {
        if (pTreeComps->GetTermNodes()[iNode] != NULL)
        {
        	pTreeComps->GetTermNodes()[iNode]->dPrediction =
                vecdNum[iNode];
            if (vecdDenom[iNode] <= 0.0)
            {
            	pTreeComps->GetTermNodes()[iNode]->dPrediction = 0.0;
            }
            else
            {
            	pTreeComps->GetTermNodes()[iNode]->dPrediction =
                    vecdNum[iNode]/vecdDenom[iNode];
            }
        }
    }
}


double CPairwise::BagImprovement
(
	const CDataset& data,
    const double *adF,
    const bag& afInBag,
    const double shrinkage,
    const double* adFadj
)
{

#ifdef NOISY_DEBUG
    Rprintf("BagImprovement, nTrain = %u\n", nTrain);
#endif

    if (data.get_trainSize() <= 0)
    {
        return 0;
    }

    double dL = 0.0;
    double dW = 0.0;

    unsigned int iItemStart = 0;
    unsigned int iItemEnd   = 0;


    while (iItemStart < data.get_trainSize())
    {
        const double dGroup = adGroup[iItemStart];

        // Find end of current group
        for (iItemEnd = iItemStart + 1; iItemEnd < data.get_trainSize() && adGroup[iItemEnd] == dGroup; iItemEnd++) ;

        if (!data.GetBagElem(iItemStart))
        {
            // Group was held out of training set

            const unsigned int cNumItems = iItemEnd - iItemStart;

            const double dMaxScore = pirm->MaxMeasure((int)dGroup, data.y_ptr() + iItemStart, cNumItems);

            if (dMaxScore > 0.0)
            {
                // If offset given, add up current scores
                const double* adFPlusOffset = OffsetVector(adF, data.offset_ptr(false), iItemStart, iItemEnd, vecdFPlusOffset);

                // Compute score according to old score, adF
                ranker.SetGroupScores(adFPlusOffset, cNumItems);
                ranker.Rank();
                const double dOldScore = pirm->Measure(data.y_ptr() + iItemStart, ranker);

                // Compute score according to new score: adF' =  adF + dStepSize * adFadj
                for (unsigned int i = 0; i < cNumItems; i++)
                {
                    ranker.AddToScore(i, adFadj[i+iItemStart] * shrinkage);
                }

                const double dWi = data.weight_ptr()[iItemStart];

                if (ranker.Rank())
                {
                    // Ranking changed
                    const double dNewScore = pirm->Measure(data.y_ptr() + iItemStart, ranker);
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
