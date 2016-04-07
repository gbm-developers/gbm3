//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//
//  File:       node.h
//
//  License:    GNU GPL (version 2 or later)
//
//  Contents:   a node in the tree
//
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//			   16/03/2016   James Hickey: updated to remove terminal and non-terminal nodes
//			   30/03/2016   James Hickey: templated to deal with continuous and categorical splits.
//------------------------------------------------------------------------------

#ifndef __node_h__
#define __node_h__
//------------------------------
// Includes
//------------------------------
#include <vector>
#include "dataset.h"
#include "nodeParameters.h"
#include "buildinfo.h"

//------------------------------
// Class Forwards and Enums
//------------------------------
class GenericNodeStrategy;
enum SplitType {categorical, continuous, none};

using namespace std;
typedef vector<int> VEC_CATEGORIES;
typedef vector<VEC_CATEGORIES> VEC_VEC_CATEGORIES;

//------------------------------
// Class definition
//------------------------------
class CNode
{
public:
	//----------------------
	// Public Constructors
	//----------------------
    CNode(double nodePrediction,
    		double trainingWeight, long numObs);

	//---------------------
	// Public destructor
	//---------------------
    virtual ~CNode();

	//---------------------
	// Public Functions
	//---------------------
    void Adjust(unsigned long cMinObsInNode);
    void Predict(const CDataset &data,
			 unsigned long iRow,
			 double &dFadj);

    void GetVarRelativeInfluence(double *adRelInf);
    void SplitNode();
    void PrintSubtree(unsigned long cIndent);
    double SplitImprovement(){ return childrenParams.ImprovedResiduals;}
    void TransferTreeToRList(int &iNodeID,
				     const CDataset &data,
				     int *aiSplitVar,
				     double *adSplitPoint,
				     int *aiLeftNode,
				     int *aiRightNode,
				     int *aiMissingNode,
				     double *adErrorReduction,
				     double *adWeight,
				     double *adPred,
				     VEC_VEC_CATEGORIES &vecSplitCodes,
				     int cCatSplitsOld,
				     double dShrinkage);
	signed char WhichNode(const CDataset &data,
							unsigned long iObs);

	//---------------------
	// Public Variables
	//---------------------
	// Pointers to the Node's children
	CNode* pLeftNode;
	CNode* pRightNode;
	CNode* pMissingNode;

	// Parameters
	NodeParams nodeParams;
	NodeParams childrenParams;

	// This nodes parameters
	unsigned long iSplitVar;
	double dImprovement;

	// Properties defining the node
	double dPrediction;
	double dTrainW;   // total training weight in node
	long cN; // number of training observations in node

	// ENUM FOR strategy
	SplitType splitType;

	// VARIABLES USED IN NODE SPLITTING
	std::vector<unsigned long> aiLeftCategory;
    double dSplitValue;

private:
	//---------------------
	// Private Functions
	//---------------------
    void SetStrategy();

	//---------------------
	// Private Variables
	//---------------------
    GenericNodeStrategy* nodeStrategy;

};

#endif // __node_h__



