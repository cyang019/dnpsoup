#ifndef DNPSOUP_PROPERTYVALUE_H
#define DNPSOUP_PROPERTYVALUE_H

#include "dnpsoup_core/spin_physics_components/spin.h"
#include <string>
#include <vector>
#include <utility>  // pair
#include <unordered_map>
#include <iostream>

namespace dnpsoup {
  enum class PropertyType : int {
    Csa = 1,
    ScalarCoupling = 2,
    DipolarCoupling = 3,
    HyperfineCoupling = 4,
    PseudoHyperfineCoupling = 5,
  };

  enum class ValueName : int {
    gyro = 0,
    bz = 1,
    sxx = 2,
    syy = 3,
    szz = 4,
    dipole = 5,
    scalar = 6,
    hyperfine = 7,
    pseudohyperfine = 8
  };

  std::ostream& operator<<(std::ostream &, const ValueName &);

  class PropertyValue {
  public:
    PropertyValue(const PropertyType &); 
    ~PropertyValue() {}

    double get(const ValueName &name) const;
    PropertyValue& set(const ValueName &name, double value);
  private:
    std::unordered_map<std::string, double> m_values;
  }; // class PropertyValue
} // namespace dnpsoup

#endif
