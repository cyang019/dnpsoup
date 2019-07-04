#ifndef DNPSOUP_ERRORS_H
#define DNPSOUP_ERRORS_H

#include <stdexcept>

namespace dnpsoup {
  class NotImplementedError : public std::logic_error {
    public:
      using std::logic_error::logic_error;
  };
} // namespace dnpsoup

#endif
