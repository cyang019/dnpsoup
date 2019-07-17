#ifndef DNPSOUP_PROPERTYVALUE_H
#define DNPSOUP_PROPERTYVALUE_H

#include "dnpsoup_core/spin_physics_components/spin.h"
#include <string>
#include <vector>
#include <utility>  // pair
#include <unordered_map>
#include <iostream>

namespace dnpsoup {
  enum class ValueName : int {
    iso = 0,  // for scalar
    xx = 10,  // for anisotropy
    yy = 11,
    zz = 12,
    xx1 = 13,  // for anisotropy
    yy1 = 14,
    zz1 = 15,
    xx2 = 16,  // for anisotropy
    yy2 = 17,
    zz2 = 18,
    offset = 100,   // for lab frame Ix Iy
    offset1 = 101,
    offset2 = 102
  };

  std::ostream& operator<<(std::ostream &, const ValueName &);

  // values are in Hz
  class PropertyValue {
  public:
    PropertyValue(const PropertyType &); 
    ~PropertyValue() {}

    double get(const ValueName &name) const;
    PropertyValue& set(const ValueName &name, double value);
  private:
    std::unordered_map<ValueName, double> m_values;
  }; // class PropertyValue

  /// @param gyro: gyromagnetic ratio in Hz/T
  /// @param bz: magnetic field (e.g. 9.4 T)
  /// @param sxx, syy, szz: chemical shift anisotropy tensor
  PropertyValue genCsaValue(double gyro, double bz, double sxx, double syy, double szz);

  /// @param gyro1, gyro2: gyromagnetic ratio of the two interacting spins.
  /// @param distance: in Anstrom
  PropertyValue genDipoleValue(double gyro1, double gyro2, double distance);

  /// @param val: value of the interaction in Hz.
  PropertyValue genScalarValue(double val);
} // namespace dnpsoup

#endif
