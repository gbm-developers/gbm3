//---------------------------------------------------------------------------------
//  GBM alteration by Stefan Schroedl (schroedl@a9.com)
//
//  File:       pairwise
//
//  Contains:   Distribution object to implement pairwise distributions for ranking
//
//  History:    12/15/2011   Created
//
//---------------------------------------------------------------------------------

//  This file implements the LambdaMart algorithm for learning ranking functions.
//  The main idea is to model p_ij, the probability that item i should rank higher
//  than j, as
//       p_ij = 1 / (1 + std::exp(s_i - s_j)),
//  where s_i, s_j are the model scores for the two items.
//
//  While scores are still generated one item at a time, gradients for learning
//  depend on _pairs_ of items. The algorithm is aware of _groups_; all pairs of items
//  with different labels, belonging to the same group, are used for training. A
//  typical application is ranking for web search: groups correspond to user queries,
//  and items to (feature vectors of) web pages in the associated match set.
//
//  Different IR measures can be chosen, to weight instances based on their rank.
//  Generally, changes in top ranks should have more influence than changes at the
//  bottom of the result list. This function provides the following options:
//
//  * CONC (concordance index, fraction of correctly raked pairs. This is a generalization
//    of Area under the ROC Curve (AUC) from binary to multivalued labels.
//  * Normalized Discounted Cumulative Gain (NDCG)
//  * Mean Reciprocal Rank (MRR) of the highest-ranked positive instance.
//  * Mean Average Precision (MAP), a generalization of MRR to multiple positive instances.
//
//  While MRR and MAP expect binary target labels, CONC and NDCG can equally work with
//  continuous values. More precisely, NDCG is defined as
//     \Sum_{r=1..n} val_r / log2(r+1),
//  where val_r is the user-specified target for the item at rank r. Note that this is
//  contrast to some definitions of NDCG that assume integer targets s_i, and
//  implicitly transform val_r = 2^{s+i}-1.
//
//  Groups are specified using an integer vector of the same length as the training instances.
//
//  Optionally, item weights can be supplied; it is assumed that all instances belonging
//  to the same group have the same weight.
//
//  For background information on LambdaMart, please see e.g. the following papers:
//
//  * Burges, C., "From RankNet to LambdaRank to LambdaMART: An Overview", Microsoft
//    Research Technical Report MSR-TR-2010-82, 2010
//  * Donmez, P., K. Svore, K.,  and Burges, C., "On the Local Optimality of
//    LambdaRank", SIGIR 2009
//  * Burges, C., Ragno, R., and Le, Q., "Learning to Rank with Non-Smooth Cost
//    Functions", NIPS 2006

#ifndef PAIRWISE_H
#define PAIRWISE_H

#include <memory>
#include "distribution.h"
#include "dataset.h"
#include "buildinfo.h"

// A class to rerank groups based on (intermediate) scores
// Note: Smaller ranks are better, the top rank is 1

class CRanker
{
public:
    // Auxiliary structure to store score and rank
    typedef std::pair<double, unsigned int> CDoubleUintPair;

    // Buffer memory allocation
    void Init(unsigned int cMaxItemsPerGroup);

    // Initialize ranker with scores of items belonging to the same group
    // - adScores is a score array, (at least) cNumItems long
    bool SetGroupScores(const double* const adScores, unsigned int cNumItems);

    // Perform the ranking
    // - Return true if any item changed its rank
    bool Rank();

    // Getter / setter
    unsigned int GetNumItems() const               { return cNumItems; }
    unsigned int GetRank(int i) const              { return vecdipScoreRank[i].second; }
    unsigned int GetItem(unsigned int iRank) const { return (vecpdipScoreRank[iRank-1] - &(vecdipScoreRank[0])); }
    void SetRank(int i, unsigned int r)            { vecdipScoreRank[i].second = r; }
    void AddToScore(int i, double delta)           { vecdipScoreRank[i].first += delta; }

protected:
    // Number of items in current group
    unsigned int cNumItems;

    // Pairs of (score, rank) for current group
    vector<CDoubleUintPair> vecdipScoreRank;

    // Array of pointers to elements of vecdipScoreRank, used for sorting
    // Note: We need a separate array for sorting in order to be able to
    // quickly look up the rank for any given item.
    vector<CDoubleUintPair*> vecpdipScoreRank;
};


// Abstract base class for all IR Measures

class CIRMeasure
{
public:
    // Constructor
    CIRMeasure() : cRankCutoff(UINT_MAX) {}

    // Destructor
    virtual ~CIRMeasure() { }

    // Getter / Setter
    unsigned int GetCutoffRank() const { return cRankCutoff; }
    void SetCutoffRank(unsigned int cRankCutoff) { this->cRankCutoff = cRankCutoff; }

    // Auxiliary function for sanity check
    bool AnyPairs(const double* const adY, unsigned int cNumItems) const
    {
        return (cNumItems >= 2                    // at least two instances
                && adY[0] > 0.0                   // at least one positive example (targets are non-increasing)
                && adY[cNumItems-1] != adY[0]);   // at least two different targets
    }

    // Memory allocation
    virtual void Init(unsigned long cMaxGroup, unsigned long cNumItems, unsigned int cRankCutoff = UINT_MAX) { this->cRankCutoff = cRankCutoff; }

     // Calculate the IR measure for the group of items set in the ranker.
     // Precondition: CRanker::SetGroupScores() has been called
     // - adY are the target scores
    virtual double Measure(const double* const adY, const CRanker& ranker) = 0;

    // Calculate the maximum achievable IR measure for a given group.
    // Side effect: the ranker state might change
    // Default implementation for MRR and MAP: if any positive items exist,
    // ranking them at the top yields a perfect measure of 1.
    virtual double MaxMeasure(unsigned int iGroup, const double* const adY, unsigned int cNumItems)
    {
        return (AnyPairs(adY, cNumItems) ? 1.0 : 0.0);
    }

    // Calculate the difference in the IR measure caused by swapping the ranks of two items.
    // Assumptions:
    // * iItemBetter has a higher label than iItemWorse (i.e., adY[iItemBetter] > adY[iItemWorse]).
    // * ranker.setGroup() has been called.
    virtual double SwapCost(int iItemBetter, int iItemWorse, const double* const adY, const CRanker& ranker) const = 0;

protected:
    // Cut-off rank below which items are ignored for measure
    unsigned int cRankCutoff;
};

// Class to implement IR Measure 'CONC' (fraction of concordant pairs). For the case of binary labels, this is
// equivalent to the area under the ROC curve (AUC).

class CConc : public CIRMeasure
{
public:
    virtual ~CConc() { }

    void Init(unsigned long cMaxGroup, unsigned long cNumItems, unsigned int cRankCutoff = UINT_MAX);

    double Measure(const double* const adY, const CRanker& ranker);

    // The maximum number of correctly classified pairs is simply all pairs with different labels
    double MaxMeasure(unsigned int iGroup, const double* const adY, unsigned int cNumItems)
    {
        return PairCount(iGroup, adY, cNumItems);
    }

    // (Cached) calculation of the number of pairs with different labels
    unsigned int PairCount(unsigned int iGroup, const double* const adY, unsigned int cNumItems);

    double SwapCost(int iItemBetter, int iItemWorse, const double* const adY, const CRanker& ranker) const;

protected:
    // Calculate the number of pairs with different labels
    int  ComputePairCount(const double* const adY, unsigned int cNumItems);

    // Caches the number of pairs with different labels, for each group
    vector<int> veccPairCount;
};

// Class to implement IR Measure 'Normalized Discounted Cumulative Gain'
// Note: Labels can have any non-negative value

class CNDCG : public CIRMeasure
{
public:

    void Init(unsigned long cMaxGroup, unsigned long cNumItems, unsigned int cRankCutoff = UINT_MAX);

    // Compute DCG
    double Measure(const double* const adY, const CRanker& ranker);

    // Compute best possible DCG
    double MaxMeasure(unsigned int iGroup, const double* const adY, unsigned int cNumItems);

    double SwapCost(int iItemBetter, int iItemWorse, const double* const adY, const CRanker& ranker) const;

protected:
     // Lookup table for rank weight (w(rank) = 1/log2(1+rank))
    vector<double> vecdRankWeight;

    // Caches the maximum achievable DCG, for each group
    vector<double> vecdMaxDCG;
};

// Class to implement IR Measure 'Mean Reciprocal Rank'
// Assumption: Labels are 0 or 1

class CMRR : public CIRMeasure
{
public:
    double Measure(const double* const adY, const CRanker& ranker);

    double SwapCost(int iItemPos, int iItemNeg, const double* const adY, const CRanker& ranker) const;

};


// Class to implement IR Measure 'Mean Average Precision'
// Assumption: Labels are 0 or 1

class CMAP : public CIRMeasure
{
public:

    void Init(unsigned long cMaxGroup, unsigned long cNumItems, unsigned int cRankCutoff = UINT_MAX);

    double Measure(const double* const adY, const CRanker& ranker);

    double SwapCost(int iItemPos, int iItemNeg, const double* const adY, const CRanker& ranker) const;
protected:

    // Buffer to hold positions of positive examples
    mutable vector<int> veccRankPos;
};


// Main class for 'pairwise' distribution
// Notes and Assumptions:
// * The items are sorted such that
//   * Instances belonging to the same group occur in
//     a contiguous range
//   * Within a group, labels are non-increasing.
// * adGroup supplies the group ID (positive integer, but double
//   format for compliance with the base class interface).
// * The targets adY are non-negative values, and binary {0,1}
//   for measures MRR and MAP.
// * Higher IR measures are better.
// * Only pairs with different labels are used for training.
// * Instance weights (adWeight) are constant among groups.
// * CPairwise::Initialize() is called before any of the other
//   functions, with same values for adY, adGroup, adWeight, and
//   nTrain. Certain values have to be precomputed for
//   efficiency.

class CPairwise : public CDistribution
{
public:

	static CDistribution* Create(SEXP radMisc,
								const char* szIRMeasure,
								int& cTrain);

    virtual ~CPairwise();

    void Initialize(const CDataset* pData);
    
    void ComputeWorkingResponse(const CDataset* pData,
    			const double *adF,
				double *adZ);
    
    double Deviance(const CDataset* pData,
    			const double *adF,
                    bool isValidationSet=false);

    double InitF(const CDataset* pData);

    void FitBestConstant(const CDataset* pData,
    		const double *adF,
			 unsigned long cTermNodes,
			 double* adZ,
			 CTreeComps* pTreeComps);

    double BagImprovement(const CDataset& data,
    					  const double *adF,
    					  const bag& afInBag,
                          const double shrinkage, const double* adFadj);

protected:

    // Constructor: determine IR measure as either "conc", "map", "mrr", or "ndcg"
    CPairwise(SEXP radMisc, const char* szIRMeasure, int& cTrain);

    // Calculate and accumulate up the gradients and Hessians from all training pairs
    void ComputeLambdas(int iGroup, unsigned int cNumItems, const double* const adY, const double* const adF, const double* const adWeight, double* adZ, double* adDeriv);

    std::auto_ptr<CIRMeasure> pirm;                 // The IR measure to use
    CRanker ranker;                   // The ranker

    vector<double> vecdHessian;       // Second derivative of loss function, for each training instance; used for Newton step

    vector<double> vecdNum;           // Buffer used for numerator   in FitBestConstant(), for each node
    vector<double> vecdDenom;         // Buffer used for denominator in FitBestConstant(), for each node

    vector<double> vecdFPlusOffset;   // Temporary buffer for (adF + adOffset), if the latter is not null
    const double* adGroup;
};

#endif // PAIRWISE_H
