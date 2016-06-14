//------------------------------------------------------------------------------
//
//  File:       varsplitter.h
//
//  Description: header for class that splits a node on a particular variable.
//
//------------------------------------------------------------------------------

#ifndef VARSPLITTER_H
#define VARSPLITTER_H

//------------------------------
// Includes
//------------------------------
#include "node.h"
#include "node_parameters.h"
#include <Rcpp.h>

//------------------------------
// Class Definition
//------------------------------
class VarSplitter {
 public:
  //----------------------
  // Public Constructors
  //----------------------
  VarSplitter();
  VarSplitter(CNode& nodetosplit, unsigned long min_num_node_obs,
   		  	  unsigned long whichvar, unsigned long numvar_classes);

  //---------------------
  // Public destructor
  //---------------------
  ~VarSplitter();

  //---------------------
  // Public Functions
  //---------------------
  void SetToSplit() { issplit_ = true; };

  void IncorporateObs(double xval, double residval, double weight,
                      long monotonicity);
  NodeParams best_split() const { return bestsplit_; }
  unsigned long SetAndReturnNumGroupMeans() {
    unsigned long num_finite_means = 0;

    for (unsigned long i = 0; i < proposedsplit_.split_class(); i++) {
      groupMeanAndCat[i].second = i;

      if (group_weight_[i] != 0.0) {
        groupMeanAndCat[i].first = group_sumresid_[i] / group_weight_[i];
        num_finite_means++;
      } else {
        groupMeanAndCat[i].first = HUGE_VAL;
      }
    }

    std::sort(groupMeanAndCat.begin(),
              groupMeanAndCat.begin() + proposedsplit_.split_class());

    return num_finite_means;
  }

  void IncrementCategories(unsigned long cat, double pred_increment,
                           double trainw_increment) {
    group_sumresid_[cat] += pred_increment;
    group_weight_[cat] += trainw_increment;
    group_num_obs_[cat]++;
  }

  void UpdateLeftNodeWithCat(long cat_index) {
    proposedsplit_.UpdateLeftNode(
        group_sumresid_[groupMeanAndCat[cat_index].second],
        group_weight_[groupMeanAndCat[cat_index].second],
        group_num_obs_[groupMeanAndCat[cat_index].second]);
  }

  void EvaluateCategoricalSplit();
  void WrapUpCurrentVariable();

 private:
  //---------------------
  // Private Variables
  //---------------------
  double initial_sumresiduals_;
  double initial_totalweight_;
  unsigned long initial_numobs_;

  bool issplit_;
  unsigned long min_num_node_obs_;
  double last_xvalue_;
  NodeParams bestsplit_, proposedsplit_;
  std::vector<double> group_sumresid_;
  std::vector<double> group_weight_;
  std::vector<unsigned long> group_num_obs_;

  // Splitting arrays for Categorical variable
  std::vector<std::pair<double, int> > groupMeanAndCat;
};
#endif  // VARSPLITTER_H
