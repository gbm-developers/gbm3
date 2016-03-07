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
#include "distribution.h"
#include "tree.h"
#include "dataset.h"
#include "node_factory.h"
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
    CTreeComps(double dLambda,
    	    unsigned long cTrain,
    	    unsigned long cFeatures,
    	    double dBagFraction,
    	    unsigned long cDepth,
    	    unsigned long cMinObsInNode,
    	    int cGroups);


	//---------------------
	// Public destructor
	//---------------------
    ~CTreeComps();

    //---------------------
	// Public Functions
	//---------------------
    void TreeInitialize(CDistribution* pDist, CNodeFactory* pNodeFact);
    void AssignTermNodes();
    void BagData(bool IsPairwise, CDistribution* pDist);
    void GrowTrees(CDistribution* pDist, int& cNodes);
    void AdjustAndShrink();
    void PredictValid(CDistribution* pDist);
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
	bag GetBag();
	std::vector<unsigned long> GetNodeAssign();
	VEC_P_NODETERMINAL GetTermNodes();

	double* GetGrad();
	double* GetRespAdj();
	const double* GetRespAdj() const;
	const double  RespAdjElem(int ind);

	double GetLambda();
	unsigned long GetTrainNo();
	unsigned long GetValidNo();

	unsigned long GetDepth();
	unsigned long GetMinNodeObs();
	int GetNoGroups();

private:
	//-------------------
	// Private Variables
	//-------------------

    // these objects are for the tree growing
    // allocate them once here for all trees to use
    bag afInBag;
    std::vector<unsigned long> aiNodeAssign;
    std::vector<CNodeSearch> aNodeSearch;
    std::auto_ptr<CCARTTree> ptreeTemp;
    VEC_P_NODETERMINAL vecpTermNodes;
    std::vector<double> adZ;
    std::vector<double> adFadj;

    double dLambda;
    unsigned long cTrain;
    unsigned long cValid;
    unsigned long cFeatures;
    unsigned long cTotalInBag;
    double dBagFraction;
    unsigned long cDepth;
    unsigned long cMinObsInNode;
    int  cGroups;
};

#endif //  __gbmTreeComps_h__
