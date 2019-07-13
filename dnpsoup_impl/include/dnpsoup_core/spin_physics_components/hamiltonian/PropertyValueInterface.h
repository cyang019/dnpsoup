#ifndef DNPSOUP_PROPERTYVALUEINTERFACE_H
#define DNPSOUP_PROPERTYVALUEINTERFACE_H

#include <vector>

namespace dnpsoup {
  class PropertyValueInterface {
  public:
    PropertyValueInterface() {}
    virtual ~PropertyValueInterface() {}

    virtual std::vector<double> get() const = 0;
    virtual PropertyValueInterface& set(const std::vector<double> &) = 0;
  }; // class PropertyValueInterface
} // namespace dnpsoup

#endif
