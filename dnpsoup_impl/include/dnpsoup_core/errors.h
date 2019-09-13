#ifndef DNPSOUP_ERRORS_H
#define DNPSOUP_ERRORS_H

#include <stdexcept>

namespace dnpsoup {
  class NotImplementedError : public std::logic_error {
  public:
    using std::logic_error::logic_error;
  };

  class IndexError : public std::logic_error {
  public:
    using std::logic_error::logic_error;
  };

  class SizeMismatchError : public std::logic_error {
  public:
    using std::logic_error::logic_error;
  };

  class PropertyNameNotFound : public std::logic_error {
  public:
    using std::logic_error::logic_error;
  };

  class DuplicationError : public std::logic_error {
  public:
    using std::logic_error::logic_error;
  };

  class InteractionTypeError : public std::logic_error {
  public:
    using std::logic_error::logic_error;
  };

  class RelaxationValueError : public std::logic_error {
  public:
    using std::logic_error::logic_error;
  };

  class PulseSequenceError : public std::logic_error {
  public:
    using std::logic_error::logic_error;
  };

  class NameNotFoundInProperty : public std::logic_error {
  public:
    using std::logic_error::logic_error;
  };
} // namespace dnpsoup

#endif
