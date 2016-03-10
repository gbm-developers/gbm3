//-----------------------------------
//
// File: gbmTreeComps.cpp
//
// Description: class that implements the public methods of the
//    gbm engine tree components.
//
//-----------------------------------

//------------------------------
// Includes
//------------------------------
#include "gbmTreeComps.h"

//----------------------------------------
// Function Members - Public
//----------------------------------------
//-----------------------------------
// Function: CTreeComps
//
// Returns: none
//
// Description: Default constructor for the gbm engine tree components.
//
// Parameters: double - the shrinkage parameter
//    unsigned long - number of training examples
//    unsigned long - number of features
//    double - fraction of training examples in bag
//    unsigned long - depth of the tree to grow
//    unsigned long - minimum number of observations in node
//    int - number of groups in data
//
//-----------------------------------
CTreeComps::CTreeComps(double dLambda,
	    unsigned long cTrain,
	    unsigned long cFeatures,
	    double dBagFraction,
	    unsigned long cDepth,
	    unsigned long cMinObsInNode,
	    int cGroups)
{
	this-> dLambda = dLambda;
	this-> cTrain = cTrain;
	this-> cFeatures = cFeatures;
	this-> dBagFraction = dBagFraction;
	this-> cDepth = cDepth;
	this-> cMinObsInNode = cMinObsInNode;
	this-> cGroups = cGroups;
	ptreeTemp.reset(new CCARTTree);
	cTotalInBag = (unsigned long) (dBagFraction * cTrain);
	cValid = 0;

	// Ensure initialization makes sense
	if (cTotalInBag <= 0)
	{
	    throw GBM::invalid_argument("you have an empty bag!");
	}
	if (cTrain <= 0)
	{
	    throw GBM::invalid_argument("you've <= 0 training instances");
	}

}

//-----------------------------------
// Function: ~CTreeComps()
//
// Returns: none
//
// Description: default destructor for the gbmTreeComps
//
// Parameters: none
//
//-----------------------------------
CTreeComps::~CTreeComps()
{
}

//-----------------------------------
// Function: Initialize
//
// Returns: none
//
// Description: initializes the tree components
//
// Parameters: CDistribution* - ptr to the distribution object in GBM
//    CNodeFactory* - ptr to the tree node factory
//
//-----------------------------------
void CTreeComps::TreeInitialize(const CDataset* pData, CNodeFactory* pNodeFact)
{
	cValid = pData->nrow() - cTrain;
	if (pData->nrow() < int(cTrain))
	{
	    throw GBM::invalid_argument("your training instances don't make sense");
	}

	adZ.assign(pData->nrow(), 0);
	adFadj.assign(pData->nrow(), 0);

	ptreeTemp->Initialize(pNodeFact);

	// array for flagging those observations in the bag
	afInBag.resize(cTrain);

	// aiNodeAssign tracks to which node each training obs belongs
	aiNodeAssign.resize(cTrain);

	// NodeSearch objects help decide which nodes to split
	aNodeSearch.resize(2 * cDepth + 1);

	for(unsigned long i=0; i<2*cDepth+1; i++)
	{
	  aNodeSearch[i].Initialize(cMinObsInNode);
	}
	vecpTermNodes.resize(2*cDepth+1, NULL);
}

//-----------------------------------
// Function: AssignTermNodes
//
// Returns: none
//
// Description: assigns null to the terminal nodes
//
// Parameters: none
//
//-----------------------------------
void CTreeComps::AssignTermNodes()
{
	vecpTermNodes.assign(2*cDepth+1, NULL);
}

//-----------------------------------
// Function: BagData
//
// Returns: none
//
// Description: put data into bags.
//
// Parameters: bool - determines if distribution is pairwise
//    CDistribution ptr - pointer to the distribution + data
//
//-----------------------------------
void CTreeComps::BagData(bool IsPairwise, CDistribution* pDist)
{
	unsigned long i = 0;
	unsigned long cBagged = 0;
	// randomly assign observations to the Bag
	if (!IsPairwise)
	{

		// regular instance based training
		for(i=0; i<cTrain && (cBagged < cTotalInBag); i++)
		{
			if(unif_rand() * (cTrain-i) < cTotalInBag - cBagged)
			{
				afInBag[i] = true;
				cBagged++;
			}
			else
			{
				afInBag[i] = false;
			}
		}

		std::fill(afInBag.begin() + i, afInBag.end(), false);
	}
	else
	{
		// for pairwise training, sampling is per group
		// therefore, we will not have exactly cTotalInBag instances

		double dLastGroup = -1;
		bool fChosen = false;
		unsigned int cBaggedGroups = 0;
		unsigned int cSeenGroups   = 0;
		unsigned int cTotalGroupsInBag = (unsigned long)(dBagFraction * cGroups);
		if (cTotalGroupsInBag <= 0)
		{
			cTotalGroupsInBag = 1;
		}

		for(i=0; i< cTrain; i++)
		{
			const double dGroup = pDist->misc_ptr(true)[i];
			if(dGroup != dLastGroup)
			{
				if (cBaggedGroups >= cTotalGroupsInBag)
				{
					break;
				}

				// Group changed, make a new decision
				fChosen = (unif_rand()*(cGroups - cSeenGroups) <
			   cTotalGroupsInBag - cBaggedGroups);
				if(fChosen)
				{
					cBaggedGroups++;
				}
				dLastGroup = dGroup;
				cSeenGroups++;
			}
			if(fChosen)
			{
				afInBag[i] = true;
				cBagged++;
			}
			else
			{
				afInBag[i] = false;
			}
		}

		// the remainder is not in the bag
		std::fill(afInBag.begin() + i, afInBag.end(), false);
	}
}

//-----------------------------------
// Function: GrowTrees
//
// Returns: none
//
// Description: grows the tree
//
// Parameters: CDistribution ptr - pointer to the distribution + data
//    int& - reference to  number of nodes in tree
//
//-----------------------------------
void CTreeComps::GrowTrees(const CDataset* pData, int& cNodes)
{
	#ifdef NOISY_DEBUG
	  Rprintf("Reset tree\n");
	#endif
	  ptreeTemp->Reset();
	#ifdef NOISY_DEBUG
	  Rprintf("grow tree\n");
	#endif

	ptreeTemp->grow(&(adZ[0]),
	                *(pData),
	                pData->weight_ptr(),
	                &(adFadj[0]),
	                cTrain,
	                cFeatures,
	                cTotalInBag,
	                dLambda,
	                cDepth,
	                cMinObsInNode,
	                afInBag,
	                aiNodeAssign,
	                 &(aNodeSearch[0]),
	                vecpTermNodes);

	#ifdef NOISY_DEBUG
	  tempTree->Print();
	#endif

	  ptreeTemp->GetNodeCount(cNodes);
	#ifdef NOISY_DEBUG
	  Rprintf("get node count=%d\n",cNodes);
	#endif
}

//-----------------------------------
// Function: AdjustAndShrink
//
// Returns: none
//
// Description: adjusts the tree and shrinks.
//
// Parameters: none
//
//-----------------------------------
void CTreeComps::AdjustAndShrink()
{
	ptreeTemp->Adjust(aiNodeAssign,
	                  &(adFadj[0]),
	                  cTrain,
	                  vecpTermNodes,
	                  cMinObsInNode);
	ptreeTemp->SetShrinkage(dLambda);
	#ifdef NOISY_DEBUG
	  ptreeTemp->Print();
	#endif
}

//-----------------------------------
// Function: TransferTreeToRList
//
// Returns: none
//
// Description:
//
// Parameters:
//
//-----------------------------------
void CTreeComps::TransferTreeToRList(const CDataset &pData,
	     int *aiSplitVar,
	     double *adSplitPoint,
	     int *aiLeftNode,
	     int *aiRightNode,
	     int *aiMissingNode,
	     double *adErrorReduction,
	     double *adWeight,
	     double *adPred,
	     VEC_VEC_CATEGORIES &vecSplitCodes,
	     int cCatSplitsOld)
{
	ptreeTemp -> TransferTreeToRList(pData,
				 aiSplitVar,
				 adSplitPoint,
				 aiLeftNode,
				 aiRightNode,
				 aiMissingNode,
				 adErrorReduction,
				 adWeight,
				 adPred,
				 vecSplitCodes,
				 cCatSplitsOld,
				 dLambda);
}

//-----------------------------------
// Function: PredictValid
//
// Returns: none
//
// Description: makes predictions on validation set
//
// Parameters: CDistribution ptr - ptr to data + dist
//
//-----------------------------------
void CTreeComps::PredictValid(const CDataset* pData)
{
	ptreeTemp->PredictValid(*(pData), cValid, &(adFadj[0]));
}

//-----------------------------------
// Function: GetBag
//
// Returns: bag
//
// Description: getter for bag
//
// Parameters: none
//
//-----------------------------------
bag CTreeComps::GetBag()
{
	return afInBag;
}

//-----------------------------------
// Function: GetNodeAssign
//
// Returns: vector<unsigned long>
//
// Description: getter node assignments
//
// Parameters: none
//
//-----------------------------------
std::vector<unsigned long> CTreeComps::GetNodeAssign()
{
	return aiNodeAssign;
}

//-----------------------------------
// Function: GetTermNodes
//
// Returns: VEC_P_NODETERMINAL
//
// Description: getter for terminal nodes
//
// Parameters: none
//
//-----------------------------------
VEC_P_NODETERMINAL CTreeComps::GetTermNodes()
{
	return vecpTermNodes;
}

//-----------------------------------
// Function: GetGrad
//
// Returns: double ptr
//
// Description: getter for ptr to loss function grad
//
// Parameters: none
//
//-----------------------------------
double* CTreeComps::GetGrad()
{
	return &(adZ[0]);
}

//-----------------------------------
// Function: GetRespAdj
//
// Returns: double ptr
//
// Description: getter for ptr to response adjustment
//
// Parameters: none
//
//-----------------------------------
double* CTreeComps::GetRespAdj()
{
	return &(adFadj[0]);
}

const double* CTreeComps::GetRespAdj() const
{
	return &(adFadj[0]);
}

//-----------------------------------
// Function: RespAdjElem
//
// Returns: const double
//
// Description: get element of response adjustment
//
// Parameters: int - index of element to get
//
//-----------------------------------
const double  CTreeComps::RespAdjElem(int ind)
{
	return adFadj[ind];
}

//-----------------------------------
// Function: GetLambda
//
// Returns: double
//
// Description: get shrinkage
//
// Parameters: none
//
//-----------------------------------
double CTreeComps::GetLambda()
{
	return dLambda;
}

//-----------------------------------
// Function: GetTrainNo
//
// Returns: unsigned long
//
// Description: get number of training examples
//
// Parameters: none
//
//-----------------------------------
unsigned long CTreeComps::GetTrainNo()
{
	return cTrain;
}

//-----------------------------------
// Function: GetValidNo
//
// Returns: unsigned long
//
// Description: get number of validation examples
//
// Parameters: none
//
//-----------------------------------
unsigned long CTreeComps::GetValidNo()
{
	return cValid;
}

//-----------------------------------
// Function: GetDepth
//
// Returns: unsigned long
//
// Description: get depth of tree
//
// Parameters: none
//
//-----------------------------------
unsigned long CTreeComps::GetDepth()
{
	return cDepth;
}

//-----------------------------------
// Function: GetMinNodeObs
//
// Returns: unsigned long
//
// Description: get min no of observation in node
//
// Parameters: none
//
//-----------------------------------
unsigned long CTreeComps::GetMinNodeObs()
{
	return cMinObsInNode;
}

//-----------------------------------
// Function: GetNoGroups
//
// Returns: int
//
// Description: get no of groups in data
//
// Parameters: none
//
//-----------------------------------
int CTreeComps::GetNoGroups()
{
	return cGroups;
}




