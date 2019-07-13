#include "dnpsoup_core/spin_physics_components/hamiltonian/properties/ChemicalShiftAnisotropyValue.h"
#include "dnpsoup_core/errors.h"
#include <vector>

using namespace std;

namespace dnpsoup {
    vector<double> ChemicalShiftAnisotropyValue::get() const 
    {
      auto res = vector<double>({m_sxx, m_syy, m_szz, m_gyro, m_bz});
      return res;
    }

    PropertyValueInterface& ChemicalShiftAnisotropyValue::set(const std::vector<double> &vals) 
    {
#ifndef NDEBUG
      if(vals.size() != 5){
        throw SizeMismatchError("Need exactly 5 arguments for xx, yy, zz, gyromagnetic ratio, and bz.");
      }
#endif
      m_sxx = vals[0];
      m_syy = vals[1];
      m_szz = vals[2];
      m_gyro = vals[3];
      m_bz = vals[4];
      return *this;
    }
} // namespace dnpsoup
