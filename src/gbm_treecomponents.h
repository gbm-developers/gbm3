//------------------------------------------------------------------------------
//
//  File:       gbmTreeComps.h
//
//  Description:   Header file for a class containing the tree components used by
//    the gbm engine.
//
//------------------------------------------------------------------------------

#ifndef GBMTREECOMPS_H
#define GBMTREECOMPS_H
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
    CTreeComps(TreeParams treeconfig);


	//---------------------
	// Public destructor
	//---------------------
    ~CTreeComps();

    //---------------------
	// Public Functions
	//---------------------
    void GrowTrees(const CDataset& kData, double* residuals, const double* kDeltaEstimate);
    void AdjustAndShrink(double* kDeltaEstimate);
    void PredictValid(const CDataset& kData, double* kDeltaEstimate);
    void TransferTreeToRList(const CDataset &kData,
		     int* splitvar,
		     double* splitvalues,
		     int* leftnodes,
		     int* rightnodes,
		     int* missingnodes,
		     double* error_reduction,
		     double* weights,
		     double* predictions,
		     VEC_VEC_CATEGORIES &splitcodes_vec,
		     int prev_categorical_splits);

    // getters
	std::vector<unsigned long>& get_node_assignments()
	{
		return data_node_assignment_;
	}
	vector<CNode*>& get_terminal_nodes()
	{
		return tree_.GetTermNodes();
	}
	const double& get_shrinkage_factor() const
	{
		return tree_.GetShrinkageConst();
	}
	const unsigned long& min_num_obs_required() const
	{
		return min_num_node_obs_;
	}
	const unsigned long& size_of_tree() const
	{
		return tree_.GetNodeCount();
	}

private:
	// Private Variables
	//-------------------

    // these objects are for the tree growing
    // allocate them once here for all trees to use
    std::vector<unsigned long> data_node_assignment_;
    CNodeSearch new_node_searcher_;
    CCARTTree tree_;

    unsigned long min_num_node_obs_;
};

#endif // GBMTREECOMPS_H
