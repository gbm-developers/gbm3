#ifndef GBMEXCEPT_H
#define GBMEXCEPT_H

#include <string>
#include <stdexcept>

namespace GBM {

  class invalid_argument : public std::runtime_error {
  public:
    invalid_argument() : std::runtime_error("invalid argument") {} ;
    invalid_argument(const std::string& msg) : std::runtime_error(msg) {};
  };

  class out_of_nodes : public std::runtime_error {
  public:
    out_of_nodes() : std::runtime_error("factory empty!") {};
  };

  class failure : public std::runtime_error {
  public:
  failure() : std::runtime_error("unspecified failure") {};
  failure(const std::string& msg) : std::runtime_error(msg) {};
  };

}

#endif // GBMEXCEPT_H
