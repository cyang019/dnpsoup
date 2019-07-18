#include "dnpsoup_core/spin_physics_components/hamiltonian/Property.h"
#include "dnpsoup_core/constants.h"
#include <sstream>

using namespace std;

namespace dnpsoup {
  std::ostream& operator<<(std::ostream &os, const ValueName &vname)
  {
    switch(vname){
      case ValueName::scalar:
        os << "scalar";
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
      case ValueName::d:
        os << "d";
        break;
      default:
        break;
    }
    return os;
  }

  Property::Property()
  {
  }
  
  double Property::get(const ValueName &name) const
  {
    return m_values.at(name);
  }
  
  Property& Property::set(const ValueName &name, double value)
  {
    m_values[name] = value;
    return *this;
  }

  Property genCsaProperty(double bz, double sxx, double syy, double szz)
  {
    Property p;
    p.set(ValueName::xx, sxx);
    p.set(ValueName::yy, syy);
    p.set(ValueName::zz, szz);
    p.set(ValueName::bz, bz);
    return p;
  }

  // distances are in Anstrom
  Property genDipoleProperty(double gyro1, double gyro2, double distance)
  {
    double val = 1.0e-7 * gyro1 * gyro2 * dnpsoup::h / (distance * distance * distance) * 1.0e30;
    Property p;
    p.set(ValueName::d, val);
    return p;
  }

  Property genScalarProperty(double val)
  {
    Property p;
    p.set(ValueName::scalar, val);
    return p;
  }
} // namespace dnpsoup
