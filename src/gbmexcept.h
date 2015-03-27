#ifndef GBMEXCEPT_H
#define GBMEXCEPT_H

#include <stdexcept>

namespace GBM {

  class invalid_argument : public std::runtime_error {
  public:
    invalid_argument() : std::runtime_error("invalid argument") {} ;
  };


}

#endif
