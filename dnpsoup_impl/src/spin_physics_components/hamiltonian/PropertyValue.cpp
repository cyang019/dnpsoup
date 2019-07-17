#include "dnpsoup_core/spin_physics_components/hamiltonian/PropertyValue.h"
#include "dnpsoup_core/constants.h"
#include <sstream>

using namespace std;

namespace dnpsoup {
  std::ostream& operator<<(std::ostream &os, const ValueName &vname);
  {
    switch(vname){
      case ValueName::iso:
        os << "iso";
        break;
      case ValueName::xx:
        os << "xx";
        break;
      case ValueName::yy:
        os << "yy";
        break;
      case ValueName::zz:
        os << "zz";
        break;
      case ValueName::offset:
        os << "offset";
        break;
      default:
        break;
      return os;
    }
  }

  PropertyValue::PropertyValue(const PropertyType &pt)
  {
  }
  
  double PropertyValue::get(const ValueName &name) const
  {
    return m_values.at(name);
  }
  
  PropertyValue& PropertyValue::set(const ValueName &name, double value)
  {
    m_values[name] = value;
    return *this;
  }

  PropertyValue genCsaValue(double gyro, double bz, double sxx, double syy, double szz)
  {
    double freq = gyro * bz;
    PropertyValue p;
    p.set(ValueName::xx, freq * sxx);
    p.set(ValueName::yy, freq * syy);
    p.set(ValueName::zz, freq * szz);
    return p;
  }

  // distances are in Anstrom
  PropertyValue genDipoleValue(double gyro1, double gyro2, double distance)
  {
    double val = 1.0e-7 * gyro1 * gyro2 * dnpsoup::h / (distance * distance * distance) * 1.0e30;
    PropertyValue p;
    p.set(ValueName::iso, val);
    return p;
  }

  PropertyValue genScalarValue(double val)
  {
    PropertyValue p;
    p.set(ValueName::iso, val);
    return p;
  }
} // namespace dnpsoup
