//------------------------------------------------------------------------------
//
//  File:       nodeParameters.h
//
//  Description:  header  contains the parameters used to split a node
//
//------------------------------------------------------------------------------

#ifndef NODEPARAMETERS_H
#define NODEPARAMETERS_H

//------------------------------
// Includes
//------------------------------
#include <Rcpp.h>

//------------------------------
// Struct Definition
//------------------------------
struct NodeDef {
  //----------------------
  // Public Constructors
  //----------------------
  NodeDef() : numobs_(0), weightresid_(0), totalweight_(0){};

  NodeDef(double weightresid, double totalweight, unsigned long numobs)
      : numobs_(numobs), weightresid_(weightresid), totalweight_(totalweight){};

  //---------------------
  // Public Functions
  //---------------------
  void clear() {
    numobs_ = 0;
    weightresid_ = totalweight_ = 0;
  };

  void increment(const double kPred, const double kTrainWeight, long num) {
    weightresid_ += kPred;
    totalweight_ += kTrainWeight;
    numobs_ += num;
  };

  double prediction() const { return weightresid_ / totalweight_; };

  double unweighted_gradient(const NodeDef& kOtherNode) const {
    return weightresid_ * kOtherNode.totalweight_ -
           kOtherNode.weightresid_ * totalweight_;
  };

  double variance_reduction(const NodeDef& kOtherNode) const {
    const double predictionDiff = prediction() - kOtherNode.prediction();
    return totalweight_ * kOtherNode.totalweight_ * predictionDiff *
           predictionDiff;
  };

  bool has_min_obs(unsigned long min_num_node_obs) const {
    return (numobs_ >= min_num_node_obs);
  }

  bool has_obs() const { return numobs_; }

  double sum_weights(const NodeDef& kOtherNode) const {
    return totalweight_ + kOtherNode.totalweight_;
  }

  double sum_weights(const NodeDef& kNodeOne, const NodeDef& kNodeTwo) const {
    return totalweight_ + kNodeOne.sum_weights(kNodeTwo);
  }

  double get_totalweight() const { return totalweight_; }

  long get_num_obs() const { return numobs_; }

 private:
  unsigned long numobs_;
  double weightresid_;
  double totalweight_;
};

//------------------------------
// Class Definition
//------------------------------
class NodeParams {
 public:
  //----------------------
  // Public Constructors
  //----------------------
  NodeParams() : category_ordering_(1024){};

  //---------------------
  // Public destructor
  //---------------------
  ~NodeParams();

  //---------------------
  // Public Functions
  //---------------------
  void ResetSplitProperties(double weightedresiduals, double trainingweight,
                            unsigned long numobs,
                            unsigned long variable_classes = 1,
                            unsigned long splitvar = UINT_MAX);
  void UpdateMissingNode(double predincrement, double trainw_increment,
                         long numincrement = 1) {
    // Move data point from right node to missing
    missing_.increment(predincrement, trainw_increment, numincrement);
    right_.increment(-predincrement, -trainw_increment, -numincrement);
  }
  void UpdateLeftNode(double predincrement, double trainw_increment,
                      long numincrement = 1) {
    // Move data point from right node to left node
    left_.increment(predincrement, trainw_increment, numincrement);
    right_.increment(-predincrement, -trainw_increment, -numincrement);
  }
  inline double get_improvement() { return improvement_; };
  bool split_is_correct_monotonicity(long specify_monotone) {
    return ((specify_monotone == 0) ||
            ((specify_monotone * right_.unweighted_gradient(left_)) > 0));
  }
  void NodeGradResiduals() {
    // Only need to look at left and right
    if (!missing_.has_obs()) {
      improvement_ =
          left_.variance_reduction(right_) / left_.sum_weights(right_);
    } else {
      // Grad - left/right
      improvement_ = (left_.variance_reduction(right_) +
                      left_.variance_reduction(missing_) +
                      right_.variance_reduction(missing_)) /
                     left_.sum_weights(right_, missing_);
    }
  };

  bool has_min_num_obs(unsigned long min_num_node_obs) {
    return (left_.has_min_obs(min_num_node_obs) &&
            right_.has_min_obs(min_num_node_obs));
  }
  inline void SetBestCategory(
      std::vector<std::pair<double, int> >& groupmean_and_cat) {
    int count = 0;
    category_ordering_.resize(groupmean_and_cat.size());
    for (std::vector<std::pair<double, int> >::const_iterator it =
             groupmean_and_cat.begin();
         it != groupmean_and_cat.end(); ++it) {
      category_ordering_[count] = it->second;
      count++;
    }
  };

  bool has_missing() const { return missing_.has_obs(); };
  //---------------------
  // Public Variables
  //---------------------
  // Left Node Definition
  NodeDef left_, right_, missing_;

  // Splitting values
  double split_value_;                  // Continuous Split Value
  unsigned long split_var_;             // Which feature to split on
  unsigned long split_class_;           // Categorical Split Value
  std::vector<int> category_ordering_;  // Vector of levels ordering
  double improvement_;
};

#endif  // NODEPARAMETERS_H
