// GBM by Greg Ridgeway  Copyright (C) 2003

#include "gbm.h"
#include <memory>
#include <Rcpp.h>

extern "C" {

SEXP gbm
(
    SEXP radY,       // outcome or response
    SEXP radOffset,  // offset for f(x), NA for no offset
    SEXP radX,
    SEXP raiXOrder,
    SEXP radWeight,
    SEXP radMisc,   // other row specific data (eg failure time), NA=no Misc
    SEXP rcRows,
    SEXP rcCols,
    SEXP racVarClasses,
    SEXP ralMonotoneVar,
    SEXP rszFamily,
    SEXP rcTrees,
    SEXP rcDepth,       // interaction depth
    SEXP rcMinObsInNode,
    SEXP rcNumClasses,
    SEXP rdShrinkage,
    SEXP rdBagFraction,
    SEXP rcTrain,
    SEXP rcFeatures,
    SEXP radFOld,
    SEXP rcCatSplitsOld,
    SEXP rcTreesOld,
    SEXP rfVerbose
)
{
  BEGIN_RCPP
    
    VEC_VEC_CATEGORIES vecSplitCodes;

    int i = 0;
    int iT = 0;
    int iK = 0;
    const int cTrees = Rcpp::as<int>(rcTrees);
    const int cTrain = Rcpp::as<int>(rcTrain);
    const int cFeatures = Rcpp::as<int>(rcFeatures);
    const int cNumClasses = Rcpp::as<int>(rcNumClasses);
    const int cDepth = Rcpp::as<int>(rcDepth);
    const int cMinObsInNode = Rcpp::as<int>(rcMinObsInNode);
    const int cCatSplitsOld = Rcpp::as<int>(rcCatSplitsOld);
    const int cTreesOld = Rcpp::as<int>(rcTreesOld);
    const double dShrinkage = Rcpp::as<double>(rdShrinkage);
    const double dBagFraction = Rcpp::as<double>(rdBagFraction);
    const bool verbose = Rcpp::as<bool>(rfVerbose);
    
    Rcpp::NumericVector adY(radY);
    Rcpp::NumericVector adOffset(radOffset);
    Rcpp::NumericVector adX(radX);
    Rcpp::IntegerVector aiXOrder(raiXOrder);
    Rcpp::NumericVector adMisc(radMisc);
    Rcpp::NumericVector adFold(radFOld);
    Rcpp::NumericVector adWeight(radWeight);
    const int cRows = Rcpp::as<int>(rcRows);
    const int cCols = Rcpp::as<int>(rcCols);
    Rcpp::IntegerVector acVarClasses(racVarClasses);
    Rcpp::IntegerVector alMonotoneVar(ralMonotoneVar);
    const std::string family = Rcpp::as<std::string>(rszFamily);
    
    int cNodes = 0;

    int cGroups = -1;

    Rcpp::RNGScope scope;
    
    // set up the dataset
    std::auto_ptr<CDataset> pData(new CDataset());
    // initialize some things
    std::auto_ptr<CDistribution> pDist(gbm_setup(adY.begin(),
						 adOffset.begin(),
						 adX.begin(),
						 aiXOrder.begin(),
						 adWeight.begin(),
						 adMisc.begin(),
						 cRows,
						 cCols,
						 acVarClasses.begin(),
						 alMonotoneVar.begin(),
						 family,
						 cTrees,
						 cDepth,
						 cMinObsInNode,
						 cNumClasses,
						 dShrinkage,
						 dBagFraction,
						 cTrain,
						 cFeatures,
						 pData.get(),
						 cGroups));
    
    std::auto_ptr<CGBM> pGBM(new CGBM());
    
    // initialize the GBM
    pGBM->Initialize(pData.get(),
		     pDist.get(),
		     dShrinkage,
		     cTrain,
		     cFeatures,
		     dBagFraction,
		     cDepth,
		     cMinObsInNode,
		     cNumClasses,
		     cGroups);

    double dInitF;
    Rcpp::NumericVector adF(pData->cRows * cNumClasses);

    pDist->Initialize(pData->adY,
		      pData->adMisc,
		      pData->adOffset,
		      pData->adWeight,
		      pData->cRows);
    
    if(ISNA(adFold[0])) // check for old predictions
    {
      // set the initial value of F as a constant
      pDist->InitF(pData->adY,
		   pData->adMisc,
		   pData->adOffset,
		   pData->adWeight,
		   dInitF,
		   cTrain);

      adF.fill(dInitF);
    }
    else
      {
	if (adFold.size() != adF.size()) {
	  throw GBM::invalid_argument("old predictions are the wrong shape");
	}

	std::copy(adFold.begin(),
		  adFold.begin() + cNumClasses * pData->cRows,
		  adF.begin());
      }

    Rcpp::NumericVector adTrainError(cTrees, 0.0);
    Rcpp::NumericVector adValidError(cTrees, 0.0);
    Rcpp::NumericVector adOOBagImprove(cTrees, 0.0);
    Rcpp::GenericVector setOfTrees(cTrees * cNumClasses);

    if(verbose)
    {
       Rprintf("Iter   TrainDeviance   ValidDeviance   StepSize   Improve\n");
    }
    for(iT=0; iT<cTrees; iT++)
      {
        // Update the parameters
        pDist->UpdateParams(adF.begin(),
			    pData->adOffset,
			    pData->adWeight,
			    cTrain);

        for (iK = 0; iK < cNumClasses; iK++)
        {
          double dTrainError = 0;
	  double dValidError = 0;
	  double dOOBagImprove = 0;
          pGBM->iterate(adF.begin(),
			dTrainError,dValidError,dOOBagImprove,
			cNodes, cNumClasses, iK);
          
          // store the performance measures
          adTrainError[iT] += dTrainError;
          adValidError[iT] += dValidError;
          adOOBagImprove[iT] += dOOBagImprove;

          Rcpp::IntegerVector iSplitVar(cNodes);
          Rcpp::NumericVector dSplitPoint(cNodes);
          Rcpp::IntegerVector iLeftNode(cNodes);
          Rcpp::IntegerVector iRightNode(cNodes);
          Rcpp::IntegerVector iMissingNode(cNodes);
          Rcpp::NumericVector dErrorReduction(cNodes);
          Rcpp::NumericVector dWeight(cNodes);
          Rcpp::NumericVector dPred(cNodes);

          gbm_transfer_to_R(pGBM.get(),
			    vecSplitCodes,
			    iSplitVar.begin(),
			    dSplitPoint.begin(),
			    iLeftNode.begin(),
			    iRightNode.begin(),
			    iMissingNode.begin(),
			    dErrorReduction.begin(),
			    dWeight.begin(),
			    dPred.begin(),
			    cCatSplitsOld);

	  setOfTrees[iK + iT * cNumClasses] = 
	    Rcpp::List::create(iSplitVar,
			       dSplitPoint,
			       iLeftNode, iRightNode, iMissingNode,
			       dErrorReduction, dWeight, dPred);
        } // Close for iK

        // print the information
        if((iT <= 9) ||
           (0 == (iT+1+cTreesOld) % 20) ||
	   (iT==cTrees-1))
        {
            R_CheckUserInterrupt();
            if(verbose)
	      {
		Rprintf("%6d %13.4f %15.4f %10.4f %9.4f\n",
			iT+1+cTreesOld,
			adTrainError[iT],
			adValidError[iT],
			dShrinkage,
			adOOBagImprove[iT]);
	      }
        }
      }

    if(verbose) Rprintf("\n");

    Rcpp::GenericVector setSplitCodes(vecSplitCodes.size());

    for(i=0; i<(int)vecSplitCodes.size(); i++)
    {
      Rcpp::IntegerVector vec(vecSplitCodes[i].size());
      std::copy(vecSplitCodes[i].begin(), vecSplitCodes[i].end(), vec.begin());
      setSplitCodes[i] = vec;
    }

    return Rcpp::List::create(dInitF,
                              adF,
                              adTrainError,
                              adValidError,
                              adOOBagImprove,
                              setOfTrees,
                              setSplitCodes);
   END_RCPP
}

SEXP gbm_pred
(
   SEXP radX,         // the data matrix
   SEXP rcRows,       // number of rows
   SEXP rcCols,       // number of columns
   SEXP rcNumClasses, // number of classes
   SEXP rcTrees,      // number of trees, may be a vector
   SEXP rdInitF,      // the initial value
   SEXP rTrees,       // the list of trees
   SEXP rCSplits,     // the list of categorical splits
   SEXP raiVarType,   // indicator of continuous/nominal
   SEXP riSingleTree  // boolean whether to return only results for one tree
)
{
   BEGIN_RCPP
   int iTree = 0;
   int iObs = 0;
   const int cRows = Rcpp::as<int>(rcRows);
   const Rcpp::IntegerVector cTrees(rcTrees);
   const Rcpp::GenericVector trees(rTrees);
   const Rcpp::IntegerVector aiVarType(raiVarType);
   const Rcpp::GenericVector cSplits(rCSplits);
   const Rcpp::NumericVector adX(radX);
   const int cCols = Rcpp::as<int>(rcCols);
   const int cNumClasses = Rcpp::as<int>(rcNumClasses);
   const bool fSingleTree = Rcpp::as<bool>(riSingleTree);
   const int cPredIterations = cTrees.size();
   int iPredIteration = 0;
   int iClass = 0;

   Rcpp::NumericVector adPredF(cRows*cNumClasses*cPredIterations);

   // initialize the predicted values
   if(!fSingleTree)
   {
     std::fill(adPredF.begin(),
               adPredF.begin() + cRows * cNumClasses,
               Rcpp::as<double>(rdInitF));
   }
   else
   {
     adPredF.fill(0.0);
   }
   iTree = 0;
   for(iPredIteration=0; iPredIteration<cTrees.size(); iPredIteration++)
   {
     const int mycTrees = cTrees[iPredIteration];
     if(fSingleTree) iTree=mycTrees-1;
     if(!fSingleTree && (iPredIteration>0))
       {
         // copy over from the last rcTrees
         std::copy(adPredF.begin() + cRows * cNumClasses * (iPredIteration -1),
                   adPredF.begin() + cRows * cNumClasses * iPredIteration,
                   adPredF.begin() + cRows * cNumClasses * iPredIteration);
       }
     while(iTree<mycTrees*cNumClasses)
       {
         for (iClass = 0; iClass < cNumClasses; iClass++)
           {
             const Rcpp::GenericVector thisTree = trees[iTree];
             const Rcpp::IntegerVector iSplitVar = thisTree[0];
             const Rcpp::NumericVector dSplitCode = thisTree[1];
             const Rcpp::IntegerVector iLeftNode = thisTree[2];
             const Rcpp::IntegerVector iRightNode = thisTree[3];
             const Rcpp::IntegerVector iMissingNode = thisTree[4];
              
             for(iObs=0; iObs<cRows; iObs++)
               {
                 int iCurrentNode = 0;
                 while(iSplitVar[iCurrentNode] != -1)
                   {
                     const double dX = adX[iSplitVar[iCurrentNode]*cRows + iObs];
                     // missing?
                     if(ISNA(dX))
                       {
                         iCurrentNode = iMissingNode[iCurrentNode];
                       }
                     // continuous?
                     else if (aiVarType[iSplitVar[iCurrentNode]] == 0)
                       {
                         if(dX < dSplitCode[iCurrentNode])
                           {
                             iCurrentNode = iLeftNode[iCurrentNode];
                           }
                         else
                           {
                             iCurrentNode = iRightNode[iCurrentNode];
                           }
                       }
                     else // categorical
                       {
                         const Rcpp::IntegerVector mySplits = cSplits[dSplitCode[iCurrentNode]];
                         if (mySplits.size() < (int)dX + 1) {
                           iCurrentNode = iMissingNode[iCurrentNode];
                         } else {
                           const int iCatSplitIndicator = mySplits[(int)dX];
                           if(iCatSplitIndicator==-1)
                             {
                               iCurrentNode = iLeftNode[iCurrentNode];
                             }
                           else if (iCatSplitIndicator==1)
                             {
                               iCurrentNode = iRightNode[iCurrentNode];
                             }
                           else // categorical level not present in training
                             {
                               iCurrentNode = iMissingNode[iCurrentNode];
                             }
                         }
                       }
                   }
                 adPredF[cRows*cNumClasses*iPredIteration+cRows*iClass+iObs] += dSplitCode[iCurrentNode]; // add the prediction
               } // iObs
             iTree++;
           } // iClass
       } // iTree
   }  // iPredIteration
   
   return Rcpp::wrap(adPredF);
   END_RCPP
}


SEXP gbm_plot
(
    SEXP radX,          // vector or matrix of points to make predictions
    SEXP rcRows,        // number of rows in X
    SEXP rcCols,        // number of columns in X
    SEXP rcNumClasses,  // number of classes
    SEXP raiWhichVar,   // length=cCols, index of which var cols of X are
    SEXP rcTrees,       // number of trees to use
    SEXP rdInitF,       // initial value
    SEXP rTrees,        // tree list object
    SEXP rCSplits,      // categorical split list object
    SEXP raiVarType     // vector of variable types
)
{
    BEGIN_RCPP
    int i = 0;
    int iTree = 0;
    int iObs = 0;
    int iClass = 0;
    const int cRows = Rcpp::as<int>(rcRows);
    const int cCols = Rcpp::as<int>(rcCols);
    const int cTrees = Rcpp::as<int>(rcTrees);
    const int cNumClasses = Rcpp::as<int>(rcNumClasses);
    const Rcpp::NumericVector adX(radX);
    const Rcpp::IntegerVector aiWhichVar(raiWhichVar);
    const Rcpp::GenericVector trees(rTrees);
    const Rcpp::GenericVector cSplits(rCSplits);
    const Rcpp::IntegerVector aiVarType(raiVarType);

    int aiNodeStack[40];
    double adWeightStack[40];

    Rcpp::NumericVector adPredF(cRows * cNumClasses,
                                Rcpp::as<double>(rdInitF));
    for(iTree=0; iTree<cTrees; iTree++)
    {
        for (iClass = 0; iClass < cNumClasses; iClass++)
        {
          const Rcpp::GenericVector thisTree = trees[iClass + iTree*cNumClasses];
          const Rcpp::IntegerVector iSplitVar = thisTree[0];
          const Rcpp::NumericVector dSplitCode = thisTree[1];
          const Rcpp::IntegerVector iLeftNode = thisTree[2];
          const Rcpp::IntegerVector iRightNode = thisTree[3];
          const Rcpp::IntegerVector iMissingNode = thisTree[4];
          const Rcpp::NumericVector dW = thisTree[6];
          for(iObs=0; iObs<cRows; iObs++)
            {
              aiNodeStack[0] = 0;
              adWeightStack[0] = 1.0;
              int cStackNodes = 1;
              while(cStackNodes > 0)
                {
                  cStackNodes--;
                  int iCurrentNode = aiNodeStack[cStackNodes];
                  
                  if(iSplitVar[iCurrentNode] == -1) // terminal node
                    {
                      adPredF[iClass*cRows + iObs] +=
                         adWeightStack[cStackNodes]*dSplitCode[iCurrentNode];
                    }
                    else // non-terminal node
                    {
                        // is this a split variable that interests me?
                        int iPredVar = -1;
                        for(i=0; (iPredVar == -1) && (i < cCols); i++)
                        {
                          if(aiWhichVar[i] == iSplitVar[iCurrentNode])
                            {
                              iPredVar = i; // split is on one that interests me
                            }
                        }

                        if(iPredVar != -1) // this split is among raiWhichVar
                        {
                          const double dX = adX[iPredVar*cRows + iObs];
                          // missing?
                          if(ISNA(dX))
                            {
                              aiNodeStack[cStackNodes] = iMissingNode[iCurrentNode];
                              cStackNodes++;
                            }
                          // continuous?
                          else if(aiVarType[iSplitVar[iCurrentNode]] == 0)
                            {
                              if(dX < dSplitCode[iCurrentNode])
                                {
                                  aiNodeStack[cStackNodes] = iLeftNode[iCurrentNode];
                                  cStackNodes++;
                                }
                              else
                                {
                                  aiNodeStack[cStackNodes] = iRightNode[iCurrentNode];
                                  cStackNodes++;
                                }
                            }
                            else // categorical
                            {
                              const Rcpp::IntegerVector catSplits = cSplits[dSplitCode[iCurrentNode]];
                              
                              const int iCatSplitIndicator = catSplits[dX];
                              if(iCatSplitIndicator==-1)
                                {
                                  aiNodeStack[cStackNodes] = iLeftNode[iCurrentNode];
                                  cStackNodes++;
                                }
                              else if(iCatSplitIndicator==1)
                                {
                                  aiNodeStack[cStackNodes] = iRightNode[iCurrentNode];
                                  cStackNodes++;
                                }
                              else // handle unused level
                                {
                                  iCurrentNode = iMissingNode[iCurrentNode];
                                }
                            }
                        } // iPredVar != -1
                        else // not interested in this split, average left and right
                          {
                            aiNodeStack[cStackNodes] = iRightNode[iCurrentNode];
                            const double dCurrentW = adWeightStack[cStackNodes];
                            adWeightStack[cStackNodes] = dCurrentW *
                              dW[iRightNode[iCurrentNode]]/
                              (dW[iLeftNode[iCurrentNode]]+
                               dW[iRightNode[iCurrentNode]]);
                            cStackNodes++;
                            aiNodeStack[cStackNodes] = iLeftNode[iCurrentNode];
                            adWeightStack[cStackNodes] =
                              dCurrentW-adWeightStack[cStackNodes-1];
                            cStackNodes++;
                        }
                    } // non-terminal node
                } // while(cStackNodes > 0)
            } // iObs
        } // iClass
    } // iTree

    return Rcpp::wrap(adPredF);
    END_RCPP
} // gbm_plot

} // end extern "C"

