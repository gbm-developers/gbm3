#ifndef GBMEXCEPT_H
#define GBMEXCEPT_H

#include <string>
#include <stdexcept>

namespace gbm_exception {

class InvalidArgument : public std::runtime_error {
 public:
  InvalidArgument() : std::runtime_error("invalid argument"){};
  InvalidArgument(const std::string& msg) : std::runtime_error(msg){};
};

class OutOfNodes : public std::runtime_error {
 public:
  OutOfNodes() : std::runtime_error("factory empty!"){};
};

class Failure : public std::runtime_error {
 public:
  Failure() : std::runtime_error("unspecified failure"){};
  Failure(const std::string& msg) : std::runtime_error(msg){};
};
}

#endif  // GBMEXCEPT_H
