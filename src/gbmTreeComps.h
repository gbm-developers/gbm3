//------------------------------------------------------------------------------
//
//  File:       gbmTreeComps.h
//
//  Description:   Header file for a class containing the tree components used by
//    the gbm engine.
//
//------------------------------------------------------------------------------

#ifndef __gbmTreeComps_h__
#define __gbmTreeComps_h__
//------------------------------
// Includes
//------------------------------
#include "buildinfo.h"
#include "tree.h"
#include "dataset.h"
#include "node_search.h"
#include <vector>
#include <memory>

//------------------------------
// Class definition
//------------------------------
class CTreeComps
{
public:
	//----------------------
	// Public Constructors
	//----------------------
    CTreeComps(TreeParams treeConfig);


	//---------------------
	// Public destructor
	//---------------------
    ~CTreeComps();

    //---------------------
	// Public Functions
	//---------------------
    void GrowTrees(const CDataset& data, double* adZ, const double* adFadj);
    void AdjustAndShrink(double * adFadj);
    void PredictValid(const CDataset& data, double* adFadj);
    void TransferTreeToRList(const CDataset &pData,
		     int *aiSplitVar,
		     double *adSplitPoint,
		     int *aiLeftNode,
		     int *aiRightNode,
		     int *aiMissingNode,
		     double *adErrorReduction,
		     double *adWeight,
		     double *adPred,
		     VEC_VEC_CATEGORIES &vecSplitCodes,
		     int cCatSplitsOld);

    // getters
	std::vector<unsigned long>& GetNodeAssign()
	{
		return aiNodeAssign;
	}
	vector<CNode*>& GetTermNodes()
	{
		return tree.GetTermNodes();
	}
	const double& ShrinkageConstant() const
	{
		return tree.GetShrinkageConst();
	}
	const unsigned long& GetMinNodeObs() const
	{
		return cMinObsInNode;
	}
	const long& GetSizeOfTree() const
	{
		return tree.GetNodeCount();
	}

private:
	// Private Variables
	//-------------------

    // these objects are for the tree growing
    // allocate them once here for all trees to use
    std::vector<unsigned long> aiNodeAssign;
    CNodeSearch aNodeSearch;
    CCARTTree tree;

    unsigned long cMinObsInNode;
};

#endif //  __gbmTreeComps_h__
