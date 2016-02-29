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

  if ((cTrain <= 0) || (data.nrow() < int(cTrain))) {
    throw GBM::invalid_argument("your training instances don't make sense");
  }
  
  cTotalInBag = (unsigned long)(dBagFraction*cTrain);

  if (cTotalInBag <= 0) {
    throw GBM::invalid_argument("you have an empty bag!");
  }
  
  adZ.assign(data.nrow(), 0);
  adFadj.assign(data.nrow(), 0);
  
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
  int &cNodes
)
{
  unsigned long i = 0;
  unsigned long cBagged = 0;

  if(!fInitialized)
  {
    throw GBM::failure();
  }

  dTrainError = 0.0;
  dValidError = 0.0;
  dOOBagImprove = 0.0;

  vecpTermNodes.assign(2*cDepth+1,NULL);

  // randomly assign observations to the Bag
  if (!IsPairwise())
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
      for(i=0; i<cTrain; i++)
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


#ifdef NOISY_DEBUG
  Rprintf("Compute working response\n");
#endif

  pDist->ComputeWorkingResponse(adF,
                                &adZ[0],
                                afInBag,
                                cTrain);

#ifdef NOISY_DEBUG
  Rprintf("Reset tree\n");
#endif
  ptreeTemp->Reset();
#ifdef NOISY_DEBUG
  Rprintf("grow tree\n");
#endif

  ptreeTemp->grow(&(adZ[0]), 
                  *pData, 
                  pData->weight_ptr() ,
                  &(adFadj[0]), 
                  cTrain, 
                  cFeatures, 
                  cTotalInBag, 
                  dLambda, 
                  cDepth,
                  cMinObsInNode, 
                  afInBag, 
                  aiNodeAssign, 
                  &aNodeSearch[0],
                  vecpTermNodes);

#ifdef NOISY_DEBUG
  ptreeTemp->Print();
#endif

  ptreeTemp->GetNodeCount(cNodes);
#ifdef NOISY_DEBUG
  Rprintf("get node count=%d\n",cNodes);
#endif

  // Now I have adF, adZ, and vecpTermNodes (new node assignments)
  // Fit the best constant within each terminal node
#ifdef NOISY_DEBUG
  Rprintf("fit best constant\n");
#endif

  pDist->FitBestConstant(&adF[0],
                         &adZ[0],
                         aiNodeAssign,
                         cTrain,
                         vecpTermNodes,
                         (2*cNodes+1)/3, // number of terminal nodes
                         cMinObsInNode,
                         afInBag,
                         &adFadj[0]);

  // update training predictions
  // fill in missing nodes where N < cMinObsInNode
  ptreeTemp->Adjust(aiNodeAssign,
                    &(adFadj[0]),
                    cTrain,
                    vecpTermNodes,
                    cMinObsInNode);
  ptreeTemp->SetShrinkage(dLambda);
#ifdef NOISY_DEBUG
  ptreeTemp->Print();
#endif

  dOOBagImprove = pDist->BagImprovement(&adF[0],
                                        &adFadj[0],
                                        afInBag,
                                        dLambda,
                                        cTrain);
    
  // update the training predictions
  for(i=0; i < cTrain; i++)
  {
    adF[i] += dLambda * adFadj[i];
  }
  
  dTrainError = pDist->Deviance(adF, cTrain);

  // update the validation predictions
  ptreeTemp->PredictValid(*pData, cValid, &(adFadj[0]));

  for(i=cTrain; i < cTrain+cValid; i++)
  {
    adF[i] += adFadj[i];
  }

  dValidError = pDist->Deviance(adF + cTrain, cValid, true);
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


