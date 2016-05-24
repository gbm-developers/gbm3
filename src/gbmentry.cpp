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
   int tree_num = 0;
   int obs_num = 0;
   const Rcpp::NumericMatrix kCovarMat(covariates);
   const int kNumCovarRows = kCovarMat.nrow();
   const Rcpp::IntegerVector kTrees(num_trees);
   const Rcpp::GenericVector kFittedTrees(fitted_trees);
   const Rcpp::IntegerVector kVarType(variable_type);
   const Rcpp::GenericVector kSplits(categorical_splits);
   const bool kSingleTree = Rcpp::as<bool>(ret_single_tree_res);
   const int kPredIterations = kTrees.size();
   int prediction_iteration = 0;
   int class_num = 0;

   if ((kCovarMat.ncol() != kVarType.size())) {
     throw GBM::InvalidArgument("shape mismatch");
   }
     
   Rcpp::NumericVector predicted_func(kNumCovarRows * kPredIterations);

   // initialize the predicted values
   if(!kSingleTree)
   {
     std::fill(predicted_func.begin(),
               predicted_func.begin() + kNumCovarRows,
               Rcpp::as<double>(initial_func_est));
   }
   else
   {
     predicted_func.fill(0.0);
   }
   tree_num = 0;
   for(prediction_iteration=0; prediction_iteration<kTrees.size(); prediction_iteration++)
   {
     const int kCurrTree = kTrees[prediction_iteration];
     if(kSingleTree) tree_num=kCurrTree-1;
     if(!kSingleTree && (prediction_iteration>0))
       {
         // copy over from the last rcTrees
         std::copy(predicted_func.begin() + kNumCovarRows * (prediction_iteration -1),
                   predicted_func.begin() + kNumCovarRows * prediction_iteration,
                   predicted_func.begin() + kNumCovarRows * prediction_iteration);
       }
     while(tree_num<kCurrTree)
       {
         const Rcpp::GenericVector kThisFitTree = kFittedTrees[tree_num];
         const Rcpp::IntegerVector kThisSplitVar = kThisFitTree[0];
         const Rcpp::NumericVector kThisSplitCode = kThisFitTree[1];
         const Rcpp::IntegerVector kThisLeftNode = kThisFitTree[2];
         const Rcpp::IntegerVector kThisRightNode = kThisFitTree[3];
         const Rcpp::IntegerVector kThisMissingNode = kThisFitTree[4];
         
         for(obs_num=0; obs_num<kNumCovarRows; obs_num++)
           {
             int iCurrentNode = 0;
             while(kThisSplitVar[iCurrentNode] != -1)
               {
                 const double dX = kCovarMat[kThisSplitVar[iCurrentNode]*kNumCovarRows + obs_num];
                 // missing?
                 if(ISNA(dX))
                   {
                     iCurrentNode = kThisMissingNode[iCurrentNode];
                   }
                 // continuous?
                 else if (kVarType[kThisSplitVar[iCurrentNode]] == 0)
                   {
                     if(dX < kThisSplitCode[iCurrentNode])
                       {
                         iCurrentNode = kThisLeftNode[iCurrentNode];
                       }
                     else
                       {
                         iCurrentNode = kThisRightNode[iCurrentNode];
                       }
                   }
                 else // categorical
                   {
                     const Rcpp::IntegerVector kMySplits = kSplits[kThisSplitCode[iCurrentNode]];
                     if (kMySplits.size() < (int)dX + 1) {
                       iCurrentNode = kThisMissingNode[iCurrentNode];
                     } else {
                       const int iCatSplitIndicator = kMySplits[(int)dX];
                       if(iCatSplitIndicator==-1)
                         {
                           iCurrentNode = kThisLeftNode[iCurrentNode];
                         }
                       else if (iCatSplitIndicator==1)
                         {
                           iCurrentNode = kThisRightNode[iCurrentNode];
                         }
                       else // categorical level not present in training
                         {
                           iCurrentNode = kThisMissingNode[iCurrentNode];
                         }
                     }
                   }
               }
             predicted_func[kNumCovarRows*prediction_iteration+kNumCovarRows*class_num+obs_num] += kThisSplitCode[iCurrentNode]; // add the prediction
           } // iObs
         tree_num++;
       } // iTree
   }  // iPredIteration
   
   return Rcpp::wrap(predicted_func);
   END_RCPP
}


SEXP gbm_plot
(
    SEXP covariates,          // vector or matrix of points to make predictions
    SEXP whichvar,   // index of which var cols of X are
    SEXP num_trees,       // number of trees to use
    SEXP init_func_est,       // initial value
    SEXP fitted_trees,        // tree list object
    SEXP categorical_splits,      // categorical split list object
    SEXP var_types     // vector of variable types
)
{
    BEGIN_RCPP
    int tree_num = 0;
    int obs_num = 0;
    int class_num = 0;
    const Rcpp::NumericMatrix kCovarMat(covariates);
    const int kNumRows = kCovarMat.nrow();
    const int kNumTrees = Rcpp::as<int>(num_trees);
    const Rcpp::IntegerVector kWhichVars(whichvar);
    const Rcpp::GenericVector kFittedTrees(fitted_trees);
    const Rcpp::GenericVector kSplits(categorical_splits);
    const Rcpp::IntegerVector kVarType(var_types);

    Rcpp::NumericVector predicted_func(kNumRows,
                                Rcpp::as<double>(init_func_est));

    if (kCovarMat.ncol() != kWhichVars.size()) {
      throw GBM::InvalidArgument("shape mismatch");
    }

    for(tree_num=0; tree_num<kNumTrees; tree_num++)
    {
      const Rcpp::GenericVector kThisTree = kFittedTrees[tree_num];
      const Rcpp::IntegerVector kThisSplitVar = kThisTree[0];
      const Rcpp::NumericVector kThisSplitCode = kThisTree[1];
      const Rcpp::IntegerVector kThisLeftNode = kThisTree[2];
      const Rcpp::IntegerVector kThisRightNode = kThisTree[3];
      const Rcpp::IntegerVector kThisMissingNode = kThisTree[4];
      const Rcpp::NumericVector kThisWeight = kThisTree[6];
      for(obs_num=0; obs_num<kNumRows; obs_num++)
        {
          NodeStack stack;
          stack.push(0, 1.0);
          while( !stack.empty() )
            {
              const std::pair<int, double> top = stack.pop();
              int iCurrentNode = top.first;
              const double dWeight = top.second;
              
              if(kThisSplitVar[iCurrentNode] == -1) // terminal node
                {
                  predicted_func[class_num*kNumRows + obs_num] +=
                    dWeight * kThisSplitCode[iCurrentNode];
                }
              else // non-terminal node
                {
                  // is this a split variable that interests me?
                  const Rcpp::IntegerVector::const_iterator
                    found = std::find(kWhichVars.begin(),
                                      kWhichVars.end(),
                                      kThisSplitVar[iCurrentNode]);
                  
                  if (found != kWhichVars.end())
                    {
                      const int kPredVar = found - kWhichVars.begin();
                      const double kXValue = kCovarMat(obs_num, kPredVar);
                      // missing?
                      if(ISNA(kXValue))
                        {
                          stack.push(kThisMissingNode[iCurrentNode],dWeight);
                        }
                      // continuous?
                      else if(kVarType[kThisSplitVar[iCurrentNode]] == 0)
                        {
                          if(kXValue < kThisSplitCode[iCurrentNode])
                            {
                              stack.push(kThisLeftNode[iCurrentNode],dWeight);
                            }
                          else
                            {
                              stack.push(kThisRightNode[iCurrentNode],dWeight);
                            }
                        }
                      else // categorical
                        {
                          const Rcpp::IntegerVector kCatSplits = kSplits[kThisSplitCode[iCurrentNode]];
                          
                          const int kCatSplitIndicator = kCatSplits[kXValue];
                          if(kCatSplitIndicator==-1)
                            {
                              stack.push(kThisLeftNode[iCurrentNode],dWeight);
                            }
                          else if(kCatSplitIndicator==1)
                            {
                              stack.push(kThisRightNode[iCurrentNode],dWeight);
                            }
                          else // handle unused level
                            {
                              stack.push(kThisMissingNode[iCurrentNode],dWeight);
                            }
                        }
                    } // iPredVar != -1
                  else // not interested in this split, average left and right
                    {
                      const int kRight = kThisRightNode[iCurrentNode];
                      const int kLeft = kThisLeftNode[iCurrentNode];
                      const double kRightWeight = kThisWeight[kRight];
                      const double kLeftWeight = kThisWeight[kLeft];
                      stack.push(kRight,
                                 dWeight * kRightWeight /
                                 (kRightWeight + kLeftWeight));
                      stack.push(kLeft,
                                 dWeight * kLeftWeight /
                                 (kRightWeight + kLeftWeight));
                    }
                } // non-terminal node
            } // while(cStackNodes > 0)
        } // iObs
    } // iTree

    return Rcpp::wrap(predicted_func);
    END_RCPP
} // gbm_plot

} // end extern "C"

