#ifndef DNPSOUP_CHEMICALSHIFTANISOTROPYVALUE_H
#define DNPSOUP_CHEMICALSHIFTANISOTROPYVALUE_H

#include "dnpsoup_core/spin_physics_components/hamiltonian/PropertyValue.h"

namespace dnpsoup {
  class ChemicalShiftAnisotropyValue : public PropertyValueInterface {
  public:
    ChemicalShiftAnisotropyValue(double sxx, double syy, double szz, double gyro, double bz)
      : PropertyValueInterface(), m_sxx(sxx), m_syy(syy), m_szz(szz), m_gyro(gyro), m_bz(bz) {}
    ~ChemicalShiftAnisotropyValue() {}

    std::vector<double> get() const override;
    PropertyValueInterface& set(const std::vector<double> &) override;
  private:
    double m_sxx;
    double m_syy;
    double m_szz;

    double m_gyro;
    double m_bz;
  };
} // namespace dnpsoup

#endif
