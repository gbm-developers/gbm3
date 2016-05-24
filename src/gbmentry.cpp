// GBM by Greg Ridgeway  Copyright (C) 2003

#include "gbm_engine.h"
#include "config_structs.h"
#include <memory>
#include <utility>
#include <Rcpp.h>

namespace {
  class NodeStack {
  public:
    bool empty() const {
      return stack.empty();
    }
    
    std::pair< int, double > pop() {
      std::pair< int, double > result = stack.back();
      stack.pop_back();
      return result;
    }

    void push(int nodeIndex, double weight=0.0) {
      stack.push_back(std::make_pair(nodeIndex, weight)); 
    }
    
  private:
    std::vector< std::pair< int, double > > stack;
  };
}

extern "C" {

SEXP gbm
(
    SEXP response,       // outcome or response
    SEXP offset_vec,  // offset for f(x), NA for no offset
    SEXP covariates,
    SEXP covar_order,
    SEXP sorted_vec,
    SEXP strata_vec,
    SEXP obs_weight,
    SEXP misc,   // other row specific data (eg failure time), NA=no Misc
    SEXP prior_coeff_var, // Prior coefficient of variation in Cox PH fit.
    SEXP row_to_obs_id, // Used in the bagging process
    SEXP var_classes,
    SEXP monotonicity_vec,
    SEXP dist_family,
    SEXP num_trees,
    SEXP tree_depth,       // interaction depth
    SEXP min_num_node_obs,
    SEXP shrinkageconstant,
    SEXP fraction_inbag,
    SEXP num_rows_in_training,
    SEXP num_obs_in_training,
    SEXP number_offeatures,
    SEXP prev_func_estimate,
    SEXP prev_category_splits,
    SEXP prev_trees_fitted,
    SEXP isverbose
)
{
  BEGIN_RCPP
    
  	// Set up consts for tree fitting and transfer to R API
  	VEC_VEC_CATEGORIES splitcodes;
    const int kNumTrees = Rcpp::as<int>(num_trees);
    const int kCatSplitsOld = Rcpp::as<int>(prev_category_splits);
    const int kTreesOld = Rcpp::as<int>(prev_trees_fitted);
    const bool kIsVerbose = Rcpp::as<bool>(isverbose);
    const Rcpp::NumericVector kPrevFuncEst(prev_func_estimate);

    // Set up parameters for initialization
    ConfigStructs gbmparams (response, offset_vec, covariates,
				covar_order, sorted_vec, strata_vec, obs_weight, misc,
				prior_coeff_var, row_to_obs_id, var_classes, monotonicity_vec,
				dist_family, num_trees, tree_depth,
				min_num_node_obs, shrinkageconstant,
				fraction_inbag, num_rows_in_training, num_obs_in_training, number_offeatures);
     Rcpp::RNGScope scope;

    // Build gbm piece-by-piece
    CGBM gbm(gbmparams);

    // Set up the function estimate
    double initial_func_est = gbm.initial_function_estimate();
    Rcpp::NumericMatrix temp_covars(covariates);
    Rcpp::NumericVector func_estimate(temp_covars.nrow());


    if(ISNA(kPrevFuncEst[0])) // check for old predictions
    {
      // set the initial value of F as a constant
      func_estimate.fill(initial_func_est);
    }
    else
    {
		if (kPrevFuncEst.size() != func_estimate.size()) {
		  throw GBM::InvalidArgument("old predictions are the wrong shape");
		}
		std::copy(kPrevFuncEst.begin(), kPrevFuncEst.end(), func_estimate.begin());
     }

    Rcpp::NumericVector training_errors(kNumTrees, 0.0);
    Rcpp::NumericVector validation_errors(kNumTrees, 0.0);
    Rcpp::NumericVector outofbag_improvement(kNumTrees, 0.0);
    Rcpp::GenericVector set_of_trees(kNumTrees);

    if(kIsVerbose)
    {
       Rprintf("Iter   TrainDeviance   ValidDeviance   StepSize   Improve\n");
    }
    for(int treenum=0; treenum<kNumTrees; treenum++)
    {
    	Rcpp::checkUserInterrupt();

        double tree_training_error = 0.0;
        double tree_validation_error = 0.0;
        double tree_out_of_bagimprov = 0.0;

        gbm.FitLearner(func_estimate.begin(),
                      tree_training_error,tree_validation_error,tree_out_of_bagimprov);

        // store the performance measures
        training_errors[treenum] += tree_training_error;
        validation_errors[treenum] += tree_validation_error;
        outofbag_improvement[treenum] += tree_out_of_bagimprov;

        Rcpp::IntegerVector split_vars(gbm.size_of_fitted_tree());
        Rcpp::NumericVector split_values(gbm.size_of_fitted_tree());
        Rcpp::IntegerVector left_nodes(gbm.size_of_fitted_tree());
        Rcpp::IntegerVector right_nodes(gbm.size_of_fitted_tree());
        Rcpp::IntegerVector missing_nodes(gbm.size_of_fitted_tree());
        Rcpp::NumericVector error_reduction(gbm.size_of_fitted_tree());
        Rcpp::NumericVector weights(gbm.size_of_fitted_tree());
        Rcpp::NumericVector node_predictions(gbm.size_of_fitted_tree());


        gbm.GBMTransferTreeToRList(split_vars.begin(),
                          split_values.begin(),
                          left_nodes.begin(),
                          right_nodes.begin(),
                          missing_nodes.begin(),
                          error_reduction.begin(),
                          weights.begin(),
                          node_predictions.begin(),
                          splitcodes,
                          kCatSplitsOld);

        set_of_trees[treenum] = 
          Rcpp::List::create(split_vars,
                             split_values,
                             left_nodes, right_nodes, missing_nodes,
                             error_reduction, weights, node_predictions);

        // print the information
        if((kIsVerbose) && ((treenum <= 9) ||
			 (0 == (treenum+1+kTreesOld) % 20) ||
			 (treenum==kNumTrees-1))) {
	  Rprintf("%6d %13.4f %15.4f %10.4f %9.4f\n",
		  treenum+1+kTreesOld,
		  training_errors[treenum],
		  validation_errors[treenum],
		  gbmparams.get_tree_config().shrinkage,
		  outofbag_improvement[treenum]);
        }

      }

    if(kIsVerbose) Rprintf("\n");
    using Rcpp::_;
    return Rcpp::List::create(_["initF"]=initial_func_est,
                              _["fit"]=func_estimate,
                              _["train.error"]=training_errors,
                              _["valid.error"]=validation_errors,
                              _["oobag.improve"]=outofbag_improvement,
                              _["trees"]=set_of_trees,
                              _["c.splits"]=splitcodes);

   END_RCPP
}

SEXP gbm_pred
(
   SEXP covariates,         // the data matrix
   SEXP num_trees,      // number of trees, may be a vector
   SEXP initial_func_est,      // the initial value
   SEXP fitted_trees,       // the list of trees
   SEXP categorical_splits,     // the list of categorical splits
   SEXP variable_type,   // indicator of continuous/nominal
   SEXP ret_single_tree_res  // boolean whether to return only results for one tree
)
{
   BEGIN_RCPP
   int iTree = 0;
   int iObs = 0;
   const Rcpp::NumericMatrix kCovarMat(covariates);
   const int kNumCovarRows = kCovarMat.nrow();
   const Rcpp::IntegerVector kTrees(num_trees);
   const Rcpp::GenericVector trees(fitted_trees);
   const Rcpp::IntegerVector aiVarType(variable_type);
   const Rcpp::GenericVector cSplits(categorical_splits);
   const bool fSingleTree = Rcpp::as<bool>(ret_single_tree_res);
   const int cPredIterations = kTrees.size();
   int iPredIteration = 0;
   int iClass = 0;

   if ((kCovarMat.ncol() != aiVarType.size())) {
     throw GBM::InvalidArgument("shape mismatch");
   }
     
   Rcpp::NumericVector adPredF(kNumCovarRows * cPredIterations);

   // initialize the predicted values
   if(!fSingleTree)
   {
     std::fill(adPredF.begin(),
               adPredF.begin() + kNumCovarRows,
               Rcpp::as<double>(initial_func_est));
   }
   else
   {
     adPredF.fill(0.0);
   }
   iTree = 0;
   for(iPredIteration=0; iPredIteration<kTrees.size(); iPredIteration++)
   {
     const int mycTrees = kTrees[iPredIteration];
     if(fSingleTree) iTree=mycTrees-1;
     if(!fSingleTree && (iPredIteration>0))
       {
         // copy over from the last rcTrees
         std::copy(adPredF.begin() + kNumCovarRows * (iPredIteration -1),
                   adPredF.begin() + kNumCovarRows * iPredIteration,
                   adPredF.begin() + kNumCovarRows * iPredIteration);
       }
     while(iTree<mycTrees)
       {
         const Rcpp::GenericVector thisTree = trees[iTree];
         const Rcpp::IntegerVector iSplitVar = thisTree[0];
         const Rcpp::NumericVector dSplitCode = thisTree[1];
         const Rcpp::IntegerVector iLeftNode = thisTree[2];
         const Rcpp::IntegerVector iRightNode = thisTree[3];
         const Rcpp::IntegerVector iMissingNode = thisTree[4];
         
         for(iObs=0; iObs<kNumCovarRows; iObs++)
           {
             int iCurrentNode = 0;
             while(iSplitVar[iCurrentNode] != -1)
               {
                 const double dX = kCovarMat[iSplitVar[iCurrentNode]*kNumCovarRows + iObs];
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
             adPredF[kNumCovarRows*iPredIteration+kNumCovarRows*iClass+iObs] += dSplitCode[iCurrentNode]; // add the prediction
           } // iObs
         iTree++;
       } // iTree
   }  // iPredIteration
   
   return Rcpp::wrap(adPredF);
   END_RCPP
}


SEXP gbm_plot
(
    SEXP radX,          // vector or matrix of points to make predictions
    SEXP raiWhichVar,   // index of which var cols of X are
    SEXP rcTrees,       // number of trees to use
    SEXP rdInitF,       // initial value
    SEXP rTrees,        // tree list object
    SEXP rCSplits,      // categorical split list object
    SEXP raiVarType     // vector of variable types
)
{
    BEGIN_RCPP
    int iTree = 0;
    int iObs = 0;
    int iClass = 0;
    const Rcpp::NumericMatrix adX(radX);
    const int cRows = adX.nrow();
    const int cTrees = Rcpp::as<int>(rcTrees);
    const Rcpp::IntegerVector aiWhichVar(raiWhichVar);
    const Rcpp::GenericVector trees(rTrees);
    const Rcpp::GenericVector cSplits(rCSplits);
    const Rcpp::IntegerVector aiVarType(raiVarType);

    Rcpp::NumericVector adPredF(cRows,
                                Rcpp::as<double>(rdInitF));

    if (adX.ncol() != aiWhichVar.size()) {
      throw GBM::InvalidArgument("shape mismatch");
    }

    for(iTree=0; iTree<cTrees; iTree++)
    {
      const Rcpp::GenericVector thisTree = trees[iTree];
      const Rcpp::IntegerVector iSplitVar = thisTree[0];
      const Rcpp::NumericVector dSplitCode = thisTree[1];
      const Rcpp::IntegerVector iLeftNode = thisTree[2];
      const Rcpp::IntegerVector iRightNode = thisTree[3];
      const Rcpp::IntegerVector iMissingNode = thisTree[4];
      const Rcpp::NumericVector dW = thisTree[6];
      for(iObs=0; iObs<cRows; iObs++)
        {
          NodeStack stack;
          stack.push(0, 1.0);
          while( !stack.empty() )
            {
              const std::pair<int, double> top = stack.pop();
              int iCurrentNode = top.first;
              const double dWeight = top.second;
              
              if(iSplitVar[iCurrentNode] == -1) // terminal node
                {
                  adPredF[iClass*cRows + iObs] +=
                    dWeight * dSplitCode[iCurrentNode];
                }
              else // non-terminal node
                {
                  // is this a split variable that interests me?
                  const Rcpp::IntegerVector::const_iterator
                    found = std::find(aiWhichVar.begin(),
                                      aiWhichVar.end(),
                                      iSplitVar[iCurrentNode]);
                  
                  if (found != aiWhichVar.end())
                    {
                      const int iPredVar = found - aiWhichVar.begin();
                      const double dX = adX(iObs, iPredVar);
                      // missing?
                      if(ISNA(dX))
                        {
                          stack.push(iMissingNode[iCurrentNode],dWeight);
                        }
                      // continuous?
                      else if(aiVarType[iSplitVar[iCurrentNode]] == 0)
                        {
                          if(dX < dSplitCode[iCurrentNode])
                            {
                              stack.push(iLeftNode[iCurrentNode],dWeight);
                            }
                          else
                            {
                              stack.push(iRightNode[iCurrentNode],dWeight);
                            }
                        }
                      else // categorical
                        {
                          const Rcpp::IntegerVector catSplits = cSplits[dSplitCode[iCurrentNode]];
                          
                          const int iCatSplitIndicator = catSplits[dX];
                          if(iCatSplitIndicator==-1)
                            {
                              stack.push(iLeftNode[iCurrentNode],dWeight);
                            }
                          else if(iCatSplitIndicator==1)
                            {
                              stack.push(iRightNode[iCurrentNode],dWeight);
                            }
                          else // handle unused level
                            {
                              stack.push(iMissingNode[iCurrentNode],dWeight);
                            }
                        }
                    } // iPredVar != -1
                  else // not interested in this split, average left and right
                    {
                      const int right = iRightNode[iCurrentNode];
                      const int left = iLeftNode[iCurrentNode];
                      const double right_weight = dW[right];
                      const double left_weight = dW[left];
                      stack.push(right,
                                 dWeight * right_weight /
                                 (right_weight + left_weight));
                      stack.push(left,
                                 dWeight * left_weight /
                                 (right_weight + left_weight));
                    }
                } // non-terminal node
            } // while(cStackNodes > 0)
        } // iObs
    } // iTree

    return Rcpp::wrap(adPredF);
    END_RCPP
} // gbm_plot

} // end extern "C"

