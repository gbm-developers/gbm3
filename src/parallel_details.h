#ifndef PARALLELDETAILS_H
#define PARALLELDETAILS_H

#include "gbm_exception.h"

// simple wrapper class to conceal details of parallelization

class parallel_details {
 public:
  parallel_details() : num_threads_(1) {}
  parallel_details(int num_threads) : num_threads_(num_threads) {
    if (num_threads_ <= 0) {
      throw gbm_exception::InvalidArgument(
          "number of threads must be strictly positive");
    }
  }

  int get_num_threads() const { return num_threads_; }

 private:
  int num_threads_;
};

#endif
