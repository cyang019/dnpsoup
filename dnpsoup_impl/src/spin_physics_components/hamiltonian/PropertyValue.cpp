#include "dnpsoup_core/spin_physics_components/hamiltonian/PropertyValue.h"
#include <sstream>

using namespace std;

namespace dnpsoup {
  std::ostream& operator<<(std::ostream &os, const ValueName &vname);
  {
    switch(vname){
      case ValueName::gyro:
        os << "gyro";
        break;
      case ValueName::bz:
        os << "bz";
        break;
      case ValueName::sxx:
        os << "sxx";
        break;
      case ValueName::syy:
        os << "syy";
        break;
      case ValueName::szz:
        os << "szz";
        break;
      case ValueName::dipole:
        os << "dipole";
        break;
      case ValueName::scalar:
        os << "scalar";
        break;
      case ValueName::hyperfine:
        os << "hyperfine";
        break;
      case ValueName::pseudohyperfine:
        os << "pseudohyperfine";
        break;
      default:
        break;
      return os;
    }
  }

  PropertyValue::PropertyValue(const PropertyType &pt)
  {
    switch(pt) {
      case PropertyType::Csa:
        m_values[ValueName::sxx] = 0;
        m_values[ValueName::syy] = 0;
        m_values[ValueName::szz] = 0;
        m_values[ValueName::gyro] = 0;
        m_values[ValueName::bz] = 0;
        break;
      default:
        break;
    }
  }
  
  double PropertyValue::get(const ValueName &name) const
  {
    return m_values.at(name);
  }
  
  PropertyValue& PropertyValue::set(const ValueName &name, double value)
  {
    if(m_values.find(name) == m_values.end()){
      std::ostringstream ss_names;

      for(auto kv : m_values){
        ss_names << kv.first << " ";
      }
      string err_str = name + " not found in keys: " + ss_names.str();
      throw PropertyNameNotFound(err_str);
    } else {
      m_values[name] = value;
    }
    return *this;
  }
} // namespace dnpsoup
