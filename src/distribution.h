//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//
//  File:       distribution.h
//
//  License:    GNU GPL (version 2 or later)
//
//  Contents:   distribution object
//        	
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include "node_terminal.h"

class CDistribution
{

public:

    CDistribution();
    virtual ~CDistribution();

// In the subsequent functions, parameters have the following meaning:
// * adY      - The target
// * adMisc   - Optional auxiliary data (the precise meaning is specific to the
//              derived class)
// * adOffset - An optional offset to the score (adF)
// * adWeight - Instance training weight
// * adF      - Current score (sum of all trees generated so far)
// * adZ      - (Negative) gradient of loss function, to be predicted by tree
// * adFadj   - Output of current tree, to be added to adF
// * cLength  - Number of instances (size of vectors)
// * afInBag  - true if instance is part of training set for current tree
//              (depends on random subsampling)
// * cIdxOff  - Offset used for multi-class training (CMultinomial).

// Initialize() is called once, before training starts.
// It gives derived classes a chance for custom preparations, e.g., to allocate 
// memory or to pre-compute values that do not change between iterations.

    virtual GBMRESULT Initialize(double *adY,
                                 double *adMisc,
                                 double *adOffset,
                                 double *adWeight,
                                 unsigned long cLength) { return GBM_OK; }

// UpdateParams() is called at the start of each iteration.
// CMultinomial uses it to normalize predictions across multiple classes.

    virtual GBMRESULT UpdateParams(double *adF,
                                   double *adOffset,
                                   double *adWeight,
                                   unsigned long cLength) = 0;

// ComputeWorkingResonse() calculates the negative gradients of the
// loss function, and stores them in adZ.
    
    virtual GBMRESULT ComputeWorkingResponse(double *adY,
                                             double *adMisc,
                                             double *adOffset,
                                             double *adF,
                                             double *adZ,
                                             double *adWeight,
                                             bool *afInBag,
                                             unsigned long cLength,
	                                     int cIdxOff) = 0;

// InitF() computes the best constant prediction for all instances, and 
// stores it in dInitF.

    virtual GBMRESULT InitF(double *adY,
                            double *adMisc,
                            double *adOffset,
                            double *adWeight,
                            double &dInitF,
                            unsigned long cLength) = 0;

// Deviance() returns the value of the loss function, based on the 
// current predictions (adF).

    virtual double Deviance(double *adY,
                            double *adMisc,
                            double *adOffset,
                            double *adWeight,
                            double *adF,
                            unsigned long cLength,
	                    int cIdxOff) = 0;

// FitBestConstant() calculates and sets prediction values for all terminal nodes 
// of the tree being currently constructed.
// Assumptions: 
// * cTermNodes is the number of terminal nodes of the tree.
// * vecpTermNodes is a vector of (pointers to) the terminal nodes of the tree, of 
//   size cTermNodes.
// * aiNodeAssign is a vector of size cLength, that maps each instance to an index 
//   into vecpTermNodes for the corresponding terminal node.

    virtual GBMRESULT FitBestConstant(double *adY,
                                    double *adMisc,
                                    double *adOffset,
                                    double *adWeight,
                                    double *adF,
                                    double *adZ,
                                    unsigned long *aiNodeAssign,
                                    unsigned long cLength,
                                    VEC_P_NODETERMINAL vecpTermNodes,
                                    unsigned long cTermNodes,
                                    unsigned long cMinObsInNode,
                                    bool *afInBag,
                                    double *adFadj,
	                            int cIdxOff) = 0;

// BagImprovement() returns the incremental difference in the loss 
// function induced by scoring with (adF + dStepSize * adFAdj) instead of adF, for 
// all instances that were not part of the training set for the current tree (i.e., 
// afInBag set to false).

    virtual double BagImprovement(double *adY,
                                  double *adMisc,
                                  double *adOffset,
                                  double *adWeight,
                                  double *adF,
                                  double *adFadj,
                                  bool *afInBag,
                                  double dStepSize,
                                  unsigned long cLength) = 0;
};

typedef CDistribution *PCDistribution;
#endif // DISTRIBUTION_H



