#ifndef DNPSOUP_PROPERTY_H
#define DNPSOUP_PROPERTY_H

#include "dnpsoup_core/spin_physics_components/spin.h"
#include "dnpsoup_core/common.h"
#include <string>
#include <vector>
#include <utility>  // pair
#include <unordered_map>
#include <iostream>

namespace dnpsoup {
  enum class ValueName : int {
    scalar = 0,  // for scalar
    distance  = 1,
    offset = 2,   // for microwave frequency 
    xx = 10,  // for anisotropy
    yy = 11,
    zz = 12,
    bz = 100,
    freq = 200,
    phase = 201,
    phase0 = 202,   /// lab frame rotatin: cos(w t + phi(t)): phase0 = w * t, phase = phi(t)
  };

  std::ostream& operator<<(std::ostream &, const ValueName &);

  // values are in Hz
  class Property {
  public:
    Property();
    ~Property() {}

    double get(const ValueName &name) const;
    Property& set(const ValueName &name, double value);
  private:
    std::unordered_map<ValueName, double, HashType<ValueName>> m_values;
  }; // class PropertyValue

  /// @param gyro: gyromagnetic ratio in Hz/T
  /// @param bz: magnetic field (e.g. 9.4 T)
  /// @param sxx, syy, szz: chemical shift anisotropy tensor
  Property genCsaProperty(double bz, double sxx, double syy, double szz)

  /// @param gyro1, gyro2: gyromagnetic ratio of the two interacting spins.
  /// @param distance: in Anstrom
  Property genDipoleProperty(double distance);

  /// @param val: value of the interaction in Hz.
  Property genScalarProperty(double val);

  /// @param freq: 0.5 * gamma * B1 in Hz
  /// @param phase: phase in rad
  /// @param offset: offset from detection frequency in Hz
  Property genWaveProperty(double freq, double phase, double offset);
} // namespace dnpsoup

#endif