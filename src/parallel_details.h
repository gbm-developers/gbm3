#ifndef PARALLELDETAILS_H
#define PARALLELDETAILS_H

#include "gbm_exception.h"

// simple wrapper class to conceal details of parallelization

class parallel_details {
 public:
  parallel_details() : num_threads_(1), array_chunk_size_(1024) {}
  parallel_details(int num_threads,
		   int array_chunk_size)
    : num_threads_(num_threads),
      array_chunk_size_(array_chunk_size) {
    if (num_threads_ <= 0) {
      throw gbm_exception::InvalidArgument(
          "number of threads must be strictly positive");
    }

    if (array_chunk_size_ <= 0) {
      throw gbm_exception::InvalidArgument(
	  "array chunk size must be strictly positive");
    }
  }

  int get_num_threads() const { return num_threads_; }
  int get_array_chunk_size() const { return array_chunk_size_; }

 private:
  int num_threads_;
  int array_chunk_size_;
};

#endif
