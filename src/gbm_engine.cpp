//  GBM by Greg Ridgeway  Copyright (C) 2003
//#define NOISY_DEBUG
#include <algorithm>

#include "gbm_engine.h"

CGBM::CGBM()
{
    cDepth = 0;
    cMinObsInNode = 0;
    dBagFraction = 0.0;
    dLambda = 0.0;
    fInitialized = false;
    cTotalInBag = 0;
    cTrain = 0;
    cFeatures = 0;
    cValid = 0;

    pDist = NULL;
    pData = NULL;
}


CGBM::~CGBM()
{
}


void CGBM::Initialize
(
    const CDataset& data,
    CDistribution *pDist,
    double dLambda,
    unsigned long cTrain,
    unsigned long cFeatures,
    double dBagFraction,
    unsigned long cDepth,
    unsigned long cMinObsInNode,
    unsigned long cNumClasses,
    int cGroups
)
{
  unsigned long i=0;
  
  if(!pDist) {
    throw GBM::invalid_argument();
  }
  
  this->pData = &data;
  this->pDist = pDist;
  this->dLambda = dLambda;
  this->cTrain = cTrain;
  this->cFeatures = cFeatures;
  this->dBagFraction = dBagFraction;
  this->cDepth = cDepth;
  this->cMinObsInNode = cMinObsInNode;
  this->cGroups = cGroups;

  // allocate the tree structure
  ptreeTemp.reset(new CCARTTree);
  
  cValid = data.nrow() - cTrain;

  if ((cTrain <= 0) || (data.nrow() < cTrain)) {
    throw GBM::invalid_argument("your training instances don't make sense");
  }
  
  cTotalInBag = (unsigned long)(dBagFraction*cTrain);

  if (cTotalInBag <= 0) {
    throw GBM::invalid_argument("you have an empty bag!");
  }
  
  adZ.assign(data.nrow() * cNumClasses, 0);
  adFadj.assign(data.nrow() * cNumClasses, 0);
  
  pNodeFactory.reset(new CNodeFactory());
  pNodeFactory->Initialize(cDepth);
  ptreeTemp->Initialize(pNodeFactory.get());
  
  // array for flagging those observations in the bag
  afInBag.resize(cTrain);
  
  // aiNodeAssign tracks to which node each training obs belongs
  aiNodeAssign.resize(cTrain);
  // NodeSearch objects help decide which nodes to split
  aNodeSearch.resize(2 * cDepth + 1);
  
  for(i=0; i<2*cDepth+1; i++)
    {
      aNodeSearch[i].Initialize(cMinObsInNode);
    }
  vecpTermNodes.resize(2*cDepth+1, NULL);
  
  fInitialized = true;
}


void CGBM::iterate
(
    double *adF,
    double &dTrainError,
    double &dValidError,
    double &dOOBagImprove,
    int &cNodes,
    int cNumClasses,
    int cClassIdx
)
{
    unsigned long i = 0;
    unsigned long cBagged = 0;
    int cIdxOff = cClassIdx * (cTrain + cValid);

    if(!fInitialized)
      {
	throw GBM::failure();
      }

    dTrainError = 0.0;
    dValidError = 0.0;
    dOOBagImprove = 0.0;

    vecpTermNodes.assign(2*cDepth+1,NULL);

    // randomly assign observations to the Bag

    if (cClassIdx == 0)
    {
        if (!IsPairwise())
        {
            // regular instance based training
            for(i=0; i<cTrain; i++) /* && (cBagged < cTotalInBag); i++) */
            {
                if(unif_rand()*(cTrain-i) < cTotalInBag-cBagged)
                {
                    afInBag[i] = true;
                    cBagged++;
                }
                else
                {
                    afInBag[i] = false;
                }
/*                if (cBagged >= cTotalInBag){
                    break;
                } */
            }
	    std::fill(afInBag.begin() + i, afInBag.end(), false);
        }
        else
        {
            // for pairwise training, sampling is per group
            // therefore, we will not have exactly cTotalInBag instances

            double dLastGroup = -1;
            bool chosen = false;
            unsigned int cBaggedGroups = 0;
            unsigned int cSeenGroups   = 0;
            unsigned int cTotalGroupsInBag = (unsigned long)(dBagFraction * cGroups);
            if (cTotalGroupsInBag <= 0)
            {
                cTotalGroupsInBag = 1;
            }
            for(i=0; i<cTrain; i++)
            {
              const double dGroup = pData->misc_ptr(true)[i];
              if (dGroup != dLastGroup)
                {
                  if (cBaggedGroups >= cTotalGroupsInBag)
                    {
                      break;
                    }
                  
                  // Group changed, make a new decision
                  chosen = (unif_rand()*(cGroups - cSeenGroups) < cTotalGroupsInBag - cBaggedGroups);
                  if (chosen)
                    {
                      cBaggedGroups++;
                    }
                  dLastGroup = dGroup;
                  cSeenGroups++;
                }
              if (chosen)
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

#ifdef NOISY_DEBUG
    Rprintf("Compute working response\n");
#endif

    pDist->ComputeWorkingResponse(pData->y_ptr(),
				  pData->misc_ptr(false),
				  pData->offset_ptr(false),
				  adF,
				  &adZ[0],
				  pData->weight_ptr(),
				  afInBag,
				  cTrain,
				  cIdxOff);

#ifdef NOISY_DEBUG
    Rprintf("Reset tree\n");
#endif
    ptreeTemp->Reset();
#ifdef NOISY_DEBUG
    Rprintf("grow tree\n");
#endif

    ptreeTemp->grow(&(adZ[cIdxOff]), *pData, pData->weight_ptr() + cIdxOff,
		    &(adFadj[cIdxOff]), cTrain, cFeatures, cTotalInBag, 
		    dLambda, cDepth,
		    cMinObsInNode, afInBag, aiNodeAssign, &aNodeSearch[0],
		    vecpTermNodes);
    
    
#ifdef NOISY_DEBUG
    Rprintf("get node count\n");
#endif
    ptreeTemp->GetNodeCount(cNodes);

    // Now I have adF, adZ, and vecpTermNodes (new node assignments)
    // Fit the best constant within each terminal node
#ifdef NOISY_DEBUG
    Rprintf("fit best constant\n");
#endif

    pDist->FitBestConstant(pData->y_ptr(),
			   pData->misc_ptr(false),
			   pData->offset_ptr(false),
			   pData->weight_ptr(),
			   &adF[0],
			   &adZ[0],
			   aiNodeAssign,
			   cTrain,
			   vecpTermNodes,
			   (2*cNodes+1)/3, // number of terminal nodes
			   cMinObsInNode,
			   afInBag,
			   &adFadj[0],
			   cIdxOff);
    
    
    // update training predictions
    // fill in missing nodes where N < cMinObsInNode
    ptreeTemp->Adjust(aiNodeAssign,&(adFadj[cIdxOff]),cTrain,
		      vecpTermNodes,cMinObsInNode);
    ptreeTemp->SetShrinkage(dLambda);

    if (cClassIdx == (cNumClasses - 1))
    {
      dOOBagImprove = pDist->BagImprovement(pData->y_ptr(),
					    pData->misc_ptr(false),
					    pData->offset_ptr(false),
					    pData->weight_ptr(),
					    &adF[0],
					    &adFadj[0],
					    afInBag,
					    dLambda,
					    cTrain);
    }
    
    // update the training predictions
    for(i=0; i < cTrain; i++)
      {
        int iIdx = i + cIdxOff;
        adF[iIdx] += dLambda * adFadj[iIdx];
      }

    dTrainError = pDist->Deviance(pData->y_ptr(),
                                  pData->misc_ptr(false),
                                  pData->offset_ptr(false),
                                  pData->weight_ptr(),
                                  adF,
                                  cTrain,
                                  cIdxOff);

    // update the validation predictions
    ptreeTemp->PredictValid(*pData,cValid,&(adFadj[cIdxOff]));

    for(i=cTrain; i < cTrain+cValid; i++)
      {
        adF[i + cIdxOff] += adFadj[i + cIdxOff];
      }
    
    dValidError =
      pDist->Deviance(pData->y_ptr(),
                      pData->misc_ptr(false),
                      pData->offset_ptr(false),
                      pData->weight_ptr(),
                      adF,
                      cValid,
                      cIdxOff + cTrain);
}


void CGBM::TransferTreeToRList
(
 int *aiSplitVar,
 double *adSplitPoint,
 int *aiLeftNode,
 int *aiRightNode,
 int *aiMissingNode,
 double *adErrorReduction,
 double *adWeight,
 double *adPred,
 VEC_VEC_CATEGORIES &vecSplitCodes,
 int cCatSplitsOld
 )
{
  ptreeTemp->TransferTreeToRList(*pData,
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


