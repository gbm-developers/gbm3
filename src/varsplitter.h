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
  VarSplitter(CNode& nodetosplit, unsigned long min_num_node_obs,
	      unsigned long whichvar, unsigned long numvar_classes,
	      long monotone);
  //---------------------
  // Public Functions
  //---------------------
  void SetToSplit() { issplit_ = true; };

  void IncorporateObs(double xval, double residval, double weight);
  NodeParams best_split() const { return bestsplit_; }
  unsigned long SetAndReturnNumGroupMeans() {
    unsigned long num_finite_means = 0;

    for (unsigned long i = 0; i < proposedsplit_.split_class(); i++) {
      groupMeanAndCat[i].second = i;

      if (group_[i].get_totalweight() != 0.0) {
        groupMeanAndCat[i].first = group_[i].prediction();
        num_finite_means++;
      } else {
        groupMeanAndCat[i].first = HUGE_VAL;
      }
    }

    std::sort(groupMeanAndCat.begin(),
              groupMeanAndCat.end());

    return num_finite_means;
  }

  void IncrementCategories(unsigned long cat, double pred_increment,
                           double trainw_increment) {
    group_[cat].increment(pred_increment, trainw_increment, 1);
  }

  void UpdateLeftNodeWithCat(long cat_index) {
    const NodeDef& def = group_[groupMeanAndCat[cat_index].second];
    proposedsplit_.UpdateLeftNode(def.get_weightresid(),
				  def.get_totalweight(),
				  def.get_num_obs());
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
  long monotonicity_;
  double last_xvalue_;
  NodeParams bestsplit_, proposedsplit_;
  std::vector<NodeDef> group_;

  // Splitting arrays for Categorical variable
  std::vector<std::pair<double, int> > groupMeanAndCat;
};
#endif  // VARSPLITTER_H
